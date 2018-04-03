### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transpose
###

### S3/S4 combo for t.DelayedMatrix
t.DelayedMatrix <- function(x) aperm(x)
setMethod("t", "DelayedMatrix", t.DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by a DelayedMatrix object (typically big). Multiplication of 2
### DelayedMatrix objects is not supported.
###

.DelayedMatrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- simplify_NULL_dimnames(list(rownames(x), colnames(y)))
    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))

    ## We're going to walk along the columns so need to increase the block
    ## length so that each block is made of at least one column.
    max_block_len <- max(get_max_block_length(type(y)), nrow(y))
    spacings <- get_spacings_for_linear_capped_length_blocks(dim(y),
                                                             max_block_len)
    y_grid <- RegularArrayGrid(dim(y), spacings)
    spacings[[1L]] <- ans_dim[[1L]]
    ans_grid <- RegularArrayGrid(ans_dim, spacings)  # parallel to 'y_grid'
    nblock <- length(y_grid)  # same as 'length(ans_grid)'
    for (b in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF=FALSE)
        y_viewport <- y_grid[[b]]
        block <- as.matrix(extract_block(y, y_viewport))
        block_ans <- x %*% block
        write_block_to_sink(block_ans, sink, ans_grid[[b]])
        if (get_verbose_block_processing())
            message("OK")
    }
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

