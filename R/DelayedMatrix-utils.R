### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by a DelayedMatrix object (typically big). Multiplication of 2
### DelayedMatrix objects is not supported.
###

.BLOCK_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- simplify_NULL_dimnames(list(rownames(x), colnames(y)))
    ans_type <- typeof(vector(type(x), 1L) * vector(type(y), 1L))
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))

    y_grid <- colGrid(y)
    ans_spacings <- c(ans_dim[[1L]], ncol(y_grid[[1L]]))
    ans_grid <- RegularArrayGrid(ans_dim, ans_spacings)  # parallel to 'y_grid'
    nblock <- length(y_grid)  # same as 'length(ans_grid)'
    for (b in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF=FALSE)
        y_viewport <- y_grid[[b]]
        block <- read_block(y, y_viewport)
        block_ans <- x %*% block
        write_block(sink, ans_grid[[b]], block_ans)
        if (get_verbose_block_processing())
            message("OK")
    }
    as(sink, "DelayedArray")
}

setMethod("%*%", c("ANY", "DelayedMatrix"),
    function(x, y)
    {
        if (!is.matrix(x)) {
            if (!is.vector(x))
                stop(wmsg("matrix multiplication of a ", class(x), " object ",
                          "by a ", class(y), " object is not supported"))
            x_len <- length(x)
            y_nrow <- nrow(y)
            if (x_len != 0L && x_len != y_nrow)
                stop(wmsg("non-conformable arguments"))
            x <- matrix(x, ncol=y_nrow)
        }
        .BLOCK_mult_by_left_matrix(x, y)
    }
)

setMethod("%*%", c("DelayedMatrix", "ANY"),
    function(x, y)
    {
        if (!is.matrix(y)) {
            if (!is.vector(y))
                stop(wmsg("matrix multiplication of a ", class(x), " object ",
                          "by a ", class(y), " object is not supported"))
            y_len <- length(y)
            x_ncol <- ncol(x)
            if (y_len != 0L && y_len != x_ncol)
                stop(wmsg("non-conformable arguments"))
            y <- matrix(y, nrow=x_ncol)
        }
        t(t(y) %*% t(x))
    }
)

.BLOCK_matrix_mult <- function(x, y)
{
    stop(wmsg("Matrix multiplication of 2 DelayedMatrix derivatives is not ",
              "supported at the moment. Only matrix multiplication between ",
              "a DelayedMatrix derivative and an ordinary matrix or vector ",
              "is supported for now."))

    x_dim <- dim(x)
    y_dim <- dim(y)
    stopifnot(length(x_dim) == 2L, length(y_dim) == 2L, ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- simplify_NULL_dimnames(list(rownames(x), colnames(y)))
    ans_type <- typeof(vector(type(x), 1L) * vector(type(y), 1L))
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))
}

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), .BLOCK_matrix_mult)

setMethod("crossprod", c("DelayedMatrix", "ANY"), function(x, y) t(x) %*% y)

setMethod("crossprod", c("ANY", "DelayedMatrix"), function(x, y) t(x) %*% y)

setMethod("tcrossprod", c("DelayedMatrix", "ANY"), function(x, y) x %*% t(y))

setMethod("tcrossprod", c("ANY", "DelayedMatrix"), function(x, y) x %*% t(y))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallelized matrix multiplication
###
### We make some careful decisions about how to parallelize %*% operations.
### This is done to minimize the redundant data processing, especially when
### chunks of data can only be read in an atomic manner.
###

.define_atomic_boundaries <- function(x)
# Defines the atomic boundaries across rows or columns of 'x'.
# Atomic boundaries refer to chunks that must be read in their entirety.
# Wholly in-memory matrices have trivial boundaries, i.e., per row/column.
# Boundaries for HDF5 matrices are defined by the chunks.
{
    grid <- chunkGrid(x)
    if (is.null(grid)) {
        out <- list(row=seq_len(nrow(x)), col=seq_len(ncol(x)))
    } else {
        N <- dim(grid)
        bygrid <- dims(grid)
        out <- list(
            row=cumsum(bygrid[seq_len(N[1]),1]),
            col=cumsum(bygrid[seq_len(N[2]) * N[1],2])
        )
    }
    out
}

.split_by_dimension <- function(ndim, nworkers, boundaries) 
# Splits a dimension of the matrix into 'nworkers' jobs,
# rounded up to the nearest 'boundaries' where possible. 
{
    ideal_split <- pmin(ndim, ndim/nworkers * seq_len(nworkers))
    rounded <- findInterval(ideal_split, boundaries, left.open=TRUE) + 1L
    by_bound <- rle(rounded)
    bound_val <- by_bound$values
    bound_len <- by_bound$lengths

    # This part could be converted to C for speed, 
    # though for small numbers of workers this is probably negligible.
    limits <- integer(nworkers)
    costs <- integer(nworkers)
    current <- 0L

    for (i in seq_along(bound_val)) {
        v <- bound_val[i]
        l <- bound_len[i]

        if (l==1L) {
            cur_bound <- boundaries[v]
            limits[current + 1L] <- cur_bound
            if (current) {
                costs[current + 1L] <- cur_bound - limits[current]
            } else { 
                costs[current + 1L] <- cur_bound
            }

        } else {
            if (v==1L) {
                lower <- 0L
            } else {
                lower <- boundaries[v-1L]
            }
            cur_size <- boundaries[v] - lower
            sub_size <- ceiling(cur_size/l)
            sub_split <- pmin(cur_size, sub_size * seq_len(l))

            chosen <- current + seq_len(l)
            limits[chosen] <- sub_split + lower 
            costs[chosen] <- cur_size
        }

        current <- current + l
    }

    list(bounds=limits, costs=costs)    
}

