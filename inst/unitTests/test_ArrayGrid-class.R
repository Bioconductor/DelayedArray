#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

test_get_spacings_for_hypercube_capped_length_blocks <- function()
{
    get_spacings_for_hypercube_capped_length_blocks <-
        DelayedArray:::get_spacings_for_hypercube_capped_length_blocks

    refdim <- c(15L, 10L, 5L, 8L, 10L)

    target <- c( 3L,  3L, 3L, 3L,  3L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 243)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 323)
    checkIdentical(target, current)

    target <- c( 3L,  3L, 3L, 4L,  3L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 324)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 431)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  3L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 432)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 575)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  4L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 576)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 767)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 3L, 4L,  4L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 768)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1023)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 4L, 4L,  4L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1024)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1279)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 5L, 4L,  4L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1280)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1599)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  4L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1600)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 1999)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 2000)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 2499)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 4L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 2500)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 3124)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 5L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 3125)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 3749)
    checkIdentical(target, current)

    target <- c( 6L,  5L, 5L, 5L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 3750)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 4499)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 5L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 4500)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 5399)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  5L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 5400)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 6479)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  6L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 6480)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 7559)
    checkIdentical(target, current)

    target <- c( 7L,  6L, 5L, 6L,  6L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 7560)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 8819)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 6L,  6L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 8820)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 10289)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  6L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 10290)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 12004)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  7L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 12005)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 13719)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 8L,  7L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 13720)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 15679)
    checkIdentical(target, current)

    target <- c( 8L,  7L, 5L, 8L,  7L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 15680)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 17919)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  7L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 17920)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 20479)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  8L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 20480)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 23039)
    checkIdentical(target, current)

    target <- c( 9L,  8L, 5L, 8L,  8L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 23040)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 25919)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  8L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 25920)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 29159)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  9L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 29160)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 32399)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L,  9L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 32400)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 35999)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L, 10L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 36000)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 39999)
    checkIdentical(target, current)

    target <- c( 10L, 10L, 5L, 8L, 10L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 40000)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 43999)
    checkIdentical(target, current)

    target <- c( 11L, 10L, 5L, 8L, 10L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 44000)
    checkIdentical(target, current)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim, 47999)
    checkIdentical(target, current)

    target <- c( 14L, 10L, 5L, 8L, 10L)
    current <- get_spacings_for_hypercube_capped_length_blocks(refdim,
                                                               prod(refdim)-1)
    checkIdentical(target, current)

    current <- get_spacings_for_hypercube_capped_length_blocks(refdim,
                                                               prod(refdim))
    checkIdentical(refdim, current)
}

### We do "linear blocks" only because they are the easiest to unsplit.
.split_array_by_block <- function(x, block.maxlength)
{
    grid <- defaultGrid(x, block.maxlength,
                        chunk.grid=NULL, block.shape="linear")
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

    for (block_maxlen in c(1:7, 29:31, 39:40, 59:60, 119:120)) {
        blocks <- .split_array_by_block(a1, block_maxlen)
        current <- .unsplit_array_by_block(blocks, a1)
        checkIdentical(a1, current)

        blocks <- .split_array_by_block(A1, block_maxlen)
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

    for (block_maxlen in seq_len(length(m1) * 2L)) {
        blocks <- .split_array_by_block(m1, block_maxlen)
        current <- .unsplit_array_by_block(blocks, m1)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(M1a, block_maxlen)
        current <- .unsplit_array_by_block(blocks, M1a)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(M1b, block_maxlen)
        current <- .unsplit_array_by_block(blocks, M1b)
        checkIdentical(m1, current)

        blocks <- .split_array_by_block(tm1, block_maxlen)
        current <- .unsplit_array_by_block(blocks, tm1)
        checkIdentical(tm1, current)

        blocks <- .split_array_by_block(tM1a, block_maxlen)
        current <- .unsplit_array_by_block(blocks, tM1a)
        checkIdentical(tm1, current)

        blocks <- .split_array_by_block(tM1b, block_maxlen)
        current <- .unsplit_array_by_block(blocks, tM1b)
        checkIdentical(tm1, current)
    }
}

