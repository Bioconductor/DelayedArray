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

.define_atomic_chunks <- function(x)
# Defines atomic chunks across rows or columns of 'x'. These are subsets of 
# contiguous rows or columns that must be read in their entirety. Wholly 
# in-memory matrices have trivial boundaries, i.e., per row/column. Chunks
# for HDF5 matrices are defined by the physical layout.
#
# The end-position boundary of each chunk along the relevant dimension
# is returned, along with the chunk size (in terms of the numbers of rows
# or columns it spans) and the cost of reading the chunk. The final
# term is usually proportional to the size and can be weighted differently
# if the file is file-backed vs in memory.
{
    grid <- chunkGrid(x)
    if (is.null(grid)) {
        row_out <- list(
            bound=seq_len(nrow(x)), 
            size=rep(1L, nrow(x))
        )
        row_out$cost <- row_out$size

        col_out <- list(
            bound=seq_len(ncol(x)), 
            size=rep(1L, ncol(x))
        )
        col_out$cost <- col_out$size

    } else {
        N <- dim(grid)
        bygrid <- dims(grid)

        row_dims <- bygrid[seq_len(N[1]),1]
        row_out <- list(
            bound=cumsum(row_dims),
            size=row_dims
        )
        # file-backed penalty (relative to 1).
        # TODO: accommodate penalties for combined seeds that are only partially file-backed.
        row_out$cost <- row_out$size * 10L 

        col_dims <- bygrid[seq_len(N[2]) * N[1],2]
        col_out <- list(
            bound=cumsum(col_dims),
            size=col_dims
        )
        col_out$cost <- col_out$size * 10L
    }

    list(row=row_out, col=col_out)
}

.split_by_dimension <- function(dim_info, nworkers) 
# Splits a dimension of the matrix into 'nworkers' jobs,
# rounded up to the nearest chunk boundaries where possible. 
# Also computes the cost of the split for each job.
{
    boundaries <- dim_info$bound
    chunk_size <- dim_info$size
    chunk_cost <- dim_info$cost

    ndim <- tail(boundaries, 1) # last boundary is always the matrix end.
    ideal_split <- pmin(ndim, ndim/nworkers * seq_len(nworkers))
    rounded <- findInterval(ideal_split, boundaries, left.open=TRUE) + 1L
    by_bound <- rle(rounded)
    bound_idx <- by_bound$values
    bound_len <- by_bound$lengths

    # This part could be converted to C for speed, 
    # though for small numbers of workers this is probably negligible.
    # TODO: cost balancing across workers.
    limits <- integer(nworkers)
    costs <- integer(nworkers)
    current <- 0L

    for (i in seq_along(bound_idx)) {
        b <- bound_idx[i]
        l <- bound_len[i]

        if (l==1L) {
            # One limit assigned to one chunk; the cost is equal
            # to the sum of per-chunk costs for all chunks since
            # the previous rounded limit.
            limits[current + 1L] <- boundaries[b]
            last_b <- if (i==1L) 0L else bound_idx[i-1L]
            costs[current + 1L] <- sum(chunk_cost[seq(last_b + 1L, b)])

        } else {
            # Multiple limits assigned to one chunk; splitting
            # the chunk into subchunks for distribution.
            cur_size <- chunk_size[b]
            sub_size <- ceiling(cur_size/l)
            sub_split <- pmin(cur_size, sub_size * seq_len(l))

            chosen <- current + seq_len(l)
            limits[chosen] <- sub_split + boundaries[b] - cur_size 

            # Cost for each job is equal to the cost of the chunk. 
            costs[chosen] <- chunk_cost[b]
        }

        current <- current + l
    }

    list(split=limits, cost=costs)
}

.evaluate_split_costs <- function(split, dim_info) 
# Given a split, evaluate the costs with respect to chunk atomicity. 
# To be used when 'split' is computed from a different matrix than 'dim_info'.
{
    boundaries <- dim_info$bound
    chunk_cost <- dim_info$cost

    last_chunk <- findInterval(split, boundaries, left.open=TRUE) + 1L
    prev_chunk <- findInterval(c(1L, head(split, -1L) + 1L), boundaries, left.open=TRUE)

    cum_cost <- c(0L, cumsum(chunk_cost))
    cum_cost[last_chunk+1L] - cum_cost[prev_chunk+1L]
}

.dispatch_mult <- function(x, y, nworkers, transposed.x=FALSE, transposed.y=FALSE) 
# Decides how to split jobs across multiple workers for x %*% y.
# Assessing based on the maximum cost for a given worker under each scheme,
# with the aim being to minimize the maximum time to complete all jobs.
{
    atomic_x <- .define_atomic_chunks(x)
    nr_x <- as.double(nrow(x)) # use double to avoid overflow.
    nc_x <- as.double(ncol(x))

    atomic_y <- .define_atomic_chunks(y)
    nr_y <- as.double(nrow(y))
    nc_y <- as.double(ncol(y))

    if (transposed.x) {
        atomic_x <- list(row=atomic_x$col, col=atomic_x$row) 
        tmp <- nr_x
        nr_x <- nc_x
        nc_x <- tmp
    }

    if (transposed.y) {
        atomic_y <- list(row=atomic_y$col, col=atomic_y$row) 
        tmp <- nr_y
        nr_y <- nc_y
        nc_y <- tmp
    }

    x_by_row <- .split_by_dimension(atomic_x$row, nworkers)
    y_by_col <- .split_by_dimension(atomic_y$col, nworkers)

    # Splitting 'x' by row, and multiplying with all of 'y'.
    cost_x_by_row <- max(x_by_row$cost * nc_x) + nr_y * sum(y_by_col$cost)

    # Splitting 'y' by column and multiplying with all of 'x'.
    cost_y_by_col <- sum(x_by_row$cost) * nc_x + max(nr_y * y_by_col$cost)

    # Splitting both of them by the common dimension.
    common_x <- .split_by_dimension(atomic_x$col, nworkers)
    cost_common_x <- max(nr_x * common_x$cost + .evaluate_split_costs(common_x$split, atomic_y$row) * nc_y) 
    common_y <- .split_by_dimension(atomic_y$row, nworkers)
    cost_common_y <- max(nr_x * .evaluate_split_costs(common_y$split, atomic_x$col) + common_y$cost * nc_y)

    all_options <- c(x_by_row=cost_x_by_row, y_by_col=cost_y_by_col, common_x=cost_common_x, common_y=cost_common_y)
    choice <- names(all_options)[which.min(all_options)]
    splits <- switch(choice, x_by_row=x_by_row$split, y_by_col=y_by_col$split, common_x=common_x$split, common_y=common_y$split)
    reported <- switch(choice, x_by_row="left", y_by_col="right", "both")
    list(choice=reported, split=splits)
}

.dispatch_self <- function(x, nworkers, transposed=FALSE) 
# Decides how to split jobs across multiple workers for crossprod(x).
{
    atomic <- .define_atomic_chunks(x)
    nr <- as.double(nrow(x)) # use double to avoid overflow.
    nc <- as.double(ncol(x))

    if (transposed) {
        atomic <- list(row=atomic$col, col=atomic$row) 
        tmp <- nr
        nr <- nc
        nc <- tmp
    }

    # Splitting 'x' by column, and taking the crossproduct with all of 'x'.
    # This is the same the other way around, so we don't bother computing that.
    single <- .split_by_dimension(atomic$col, nworkers)
    cost_single <- max(nr * single$cost) + sum(single$cost) * nr

    # Splitting by the common dimension, i.e., along the rows of 'x'
    # and taking the crossproduct of each split block.
    both <- .split_by_dimension(atomic$row, nworkers)
    cost_both <- max(both$cost * nc) * 2

    all_options <- c(both=cost_both, single=cost_single)
    choice <- names(all_options)[which.min(all_options)]
    splits <- switch(choice, both=both$split, single=single$split)
    list(choice=choice, split=splits)
}
