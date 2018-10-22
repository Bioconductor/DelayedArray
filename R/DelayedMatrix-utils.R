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

