#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

test_get_max_spacings_for_hypercube_blocks <- function()
{
    get_max_spacings_for_hypercube_blocks <-
        DelayedArray:::get_max_spacings_for_hypercube_blocks

    refdim <- c(15L, 10L, 5L, 8L, 10L)

    target <- c( 3L,  3L, 3L, 3L,  3L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 243)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 323)
    checkIdentical(target, current)

    target <- c( 3L,  3L, 3L, 4L,  3L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 324)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 431)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  3L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 432)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 575)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  4L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 576)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 767)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 3L, 4L,  4L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 768)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1023)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 4L, 4L,  4L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1024)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1279)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 5L, 4L,  4L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1280)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1599)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  4L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1600)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 1999)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 2000)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 2499)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 4L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 2500)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 3124)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 5L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 3125)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 3749)
    checkIdentical(target, current)

    target <- c( 6L,  5L, 5L, 5L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 3750)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 4499)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 5L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 4500)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 5399)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  5L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 5400)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 6479)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  6L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 6480)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 7559)
    checkIdentical(target, current)

    target <- c( 7L,  6L, 5L, 6L,  6L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 7560)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 8819)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 6L,  6L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 8820)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 10289)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  6L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 10290)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 12004)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  7L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 12005)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 13719)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 8L,  7L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 13720)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 15679)
    checkIdentical(target, current)

    target <- c( 8L,  7L, 5L, 8L,  7L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 15680)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 17919)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  7L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 17920)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 20479)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  8L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 20480)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 23039)
    checkIdentical(target, current)

    target <- c( 9L,  8L, 5L, 8L,  8L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 23040)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 25919)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  8L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 25920)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 29159)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  9L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 29160)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 32399)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L,  9L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 32400)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 35999)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L, 10L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 36000)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 39999)
    checkIdentical(target, current)

    target <- c( 10L, 10L, 5L, 8L, 10L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 40000)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 43999)
    checkIdentical(target, current)

    target <- c( 11L, 10L, 5L, 8L, 10L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 44000)
    checkIdentical(target, current)
    current <- get_max_spacings_for_hypercube_blocks(refdim, 47999)
    checkIdentical(target, current)

    target <- c( 14L, 10L, 5L, 8L, 10L)
    current <- get_max_spacings_for_hypercube_blocks(refdim, prod(refdim)-1)
    checkIdentical(target, current)

    current <- get_max_spacings_for_hypercube_blocks(refdim, prod(refdim))
    checkIdentical(refdim, current)
}

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

