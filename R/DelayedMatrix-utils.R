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

### Write a new HDF5 dataset to disk. Return an HDF5Matrix object that points
### to this new dataset.
.DelayedMatrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- list(rownames(x), colnames(y))
    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))
    colblock_APPLY(y, function(submatrix) x %*% submatrix, sink=sink)
    as(sink, "DelayedArray")
}

setMethod("%*%", c("DelayedMatrix", "matrix"),
    function(x, y) t(t(y) %*% t(x))
)

setMethod("%*%", c("matrix", "DelayedMatrix"),
    .DelayedMatrix_block_mult_by_left_matrix
)

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"),
    function(x, y)
        stop(wmsg("multiplication of 2 DelayedMatrix objects is not ",
                  "supported, only multiplication of an ordinary matrix by ",
                  "a DelayedMatrix object at the moment"))
)

