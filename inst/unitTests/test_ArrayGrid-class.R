#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

test_split_array_in_linear_blocks <- function()
{
    split_array_in_linear_blocks <-
        DelayedArray:::split_array_in_linear_blocks
    unsplit_array_from_linear_blocks <-
        DelayedArray:::unsplit_array_from_linear_blocks

    a1 <- array(1:300, c(3, 10, 2, 5))
    A1 <- realize(a1)

    for (max_block_len in c(1:7, 29:31, 39:40, 59:60, 119:120)) {
        blocks <- split_array_in_linear_blocks(a1, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, a1)
        checkIdentical(a1, current)

        blocks <- split_array_in_linear_blocks(A1, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, A1)
        checkIdentical(a1, current)
    }
}

test_split_matrix_in_blocks <- function()
{
    split_array_in_linear_blocks <-
        DelayedArray:::split_array_in_linear_blocks
    unsplit_array_from_linear_blocks <-
        DelayedArray:::unsplit_array_from_linear_blocks

    a1 <- array(1:300, c(3, 10, 2, 5))
    A1 <- realize(a1)

    m1 <- a1[2, c(9, 3:7), 2, -4]
    M1a <- drop(A1[2, c(9, 3:7), 2, -4])
    checkIdentical(m1, as.matrix(M1a))

    M1b <- realize(m1)
    checkIdentical(m1, as.matrix(M1b))

    tm1 <- t(m1)
    tM1a <- t(M1a)
    checkIdentical(tm1, as.matrix(tM1a))

    tM1b <- t(M1b)
    checkIdentical(tm1, as.matrix(tM1b))

    for (max_block_len in seq_len(length(m1) * 2L)) {
        blocks <- split_array_in_linear_blocks(m1, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, m1)
        checkIdentical(m1, current)

        blocks <- split_array_in_linear_blocks(M1a, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, M1a)
        checkIdentical(m1, current)

        blocks <- split_array_in_linear_blocks(M1b, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, M1b)
        checkIdentical(m1, current)

        blocks <- split_array_in_linear_blocks(tm1, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, tm1)
        checkIdentical(tm1, current)

        blocks <- split_array_in_linear_blocks(tM1a, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, tM1a)
        checkIdentical(tm1, current)

        blocks <- split_array_in_linear_blocks(tM1b, max_block_len)
        current <- unsplit_array_from_linear_blocks(blocks, tM1b)
        checkIdentical(tm1, current)
    }
}

