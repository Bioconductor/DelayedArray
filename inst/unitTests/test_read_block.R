#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

### We do "linear blocks" (i.e. block.shape="first-dim-grows-first") only,
### because they are the easiest to unsplit.
.split_array_by_block <- function(x, block.length)
{
    grid <- defaultAutoGrid(x, block.length, chunk.grid=NULL,
                               block.shape="first-dim-grows-first")
    lapply(grid, function(viewport) read_block(x, viewport))
}

### A simple unsplit() that works only because the blocks are assumed to
### be "linear".
.unsplit_array_by_block <- function(blocks, x)
{
    ans <- DelayedArray:::combine_array_objects(blocks)
    DelayedArray:::set_dim(ans, dim(x))
}

test_split_and_unsplit_array <- function()
{
    a1 <- array(1:300, c(3, 10, 2, 5))
    A1 <- realize(a1)

    for (block_len in c(1:7, 29:31, 39:40, 59:60, 119:120)) {
        blocks <- .split_array_by_block(a1, block_len)
        current <- .unsplit_array_by_block(blocks, a1)
        checkIdentical(a1, current)

        blocks <- .split_array_by_block(A1, block_len)
        current <- .unsplit_array_by_block(blocks, A1)
        checkIdentical(a1, current)
    }
}

test_split_and_unsplit_matrix <- function()
{
    a1 <- array(1:300, c(3, 10, 2, 5))
    A1 <- realize(a1)

    m1 <- a1[2, c(9, 3:7), 2, -4]
    M1a <- A1[2, c(9, 3:7), 2, -4]
    checkIdentical(m1, as.matrix(M1a))

    M1b <- realize(m1)
    checkIdentical(m1, as.matrix(M1b))

    tm1 <- t(m1)
    tM1a <- t(M1a)
    checkIdentical(tm1, as.matrix(tM1a))

    tM1b <- t(M1b)
    checkIdentical(tm1, as.matrix(tM1b))

    for (block_len in seq_len(length(m1) * 2L)) {
        blocks <- .split_array_by_block(m1, block_len)
        current <- .unsplit_array_by_block(blocks, m1)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(M1a, block_len)
        current <- .unsplit_array_by_block(blocks, M1a)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(M1b, block_len)
        current <- .unsplit_array_by_block(blocks, M1b)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(tm1, block_len)
        current <- .unsplit_array_by_block(blocks, tm1)
        checkIdentical(tm1, current)

        blocks <- .split_array_by_block(tM1a, block_len)
        current <- .unsplit_array_by_block(blocks, tM1a)
        checkIdentical(tm1, current)

        blocks <- .split_array_by_block(tM1b, block_len)
        current <- .unsplit_array_by_block(blocks, tM1b)
        checkIdentical(tm1, current)
    }
}

