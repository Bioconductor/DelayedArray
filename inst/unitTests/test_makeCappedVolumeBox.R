#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

test_make_capped_volume_hypercube_box <- function()
{
    make_capped_volume_hypercube_box <-
        DelayedArray:::.make_capped_volume_hypercube_box

    refdim <- c(15L, 10L, 5L, 8L, 10L)

    target <- c( 3L,  3L, 3L, 3L,  3L)
    current <- make_capped_volume_hypercube_box(243L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(323L, refdim)
    checkIdentical(target, current)

    target <- c( 3L,  3L, 3L, 4L,  3L)
    current <- make_capped_volume_hypercube_box(324L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(431L, refdim)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  3L)
    current <- make_capped_volume_hypercube_box(432L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(575L, refdim)
    checkIdentical(target, current)

    target <- c( 3L,  4L, 3L, 4L,  4L)
    current <- make_capped_volume_hypercube_box(576L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(767L, refdim)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 3L, 4L,  4L)
    current <- make_capped_volume_hypercube_box(768L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(1023L, refdim)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 4L, 4L,  4L)
    current <- make_capped_volume_hypercube_box(1024L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(1279L, refdim)
    checkIdentical(target, current)

    target <- c( 4L,  4L, 5L, 4L,  4L)
    current <- make_capped_volume_hypercube_box(1280L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(1599L, refdim)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  4L)
    current <- make_capped_volume_hypercube_box(1600L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(1999L, refdim)
    checkIdentical(target, current)

    target <- c( 4L,  5L, 5L, 4L,  5L)
    current <- make_capped_volume_hypercube_box(2000L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(2499L, refdim)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 4L,  5L)
    current <- make_capped_volume_hypercube_box(2500L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(3124L, refdim)
    checkIdentical(target, current)

    target <- c( 5L,  5L, 5L, 5L,  5L)
    current <- make_capped_volume_hypercube_box(3125L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(3749L, refdim)
    checkIdentical(target, current)

    target <- c( 6L,  5L, 5L, 5L,  5L)
    current <- make_capped_volume_hypercube_box(3750L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(4499L, refdim)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 5L,  5L)
    current <- make_capped_volume_hypercube_box(4500L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(5399L, refdim)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  5L)
    current <- make_capped_volume_hypercube_box(5400L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(6479L, refdim)
    checkIdentical(target, current)

    target <- c( 6L,  6L, 5L, 6L,  6L)
    current <- make_capped_volume_hypercube_box(6480L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(7559L, refdim)
    checkIdentical(target, current)

    target <- c( 7L,  6L, 5L, 6L,  6L)
    current <- make_capped_volume_hypercube_box(7560L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(8819L, refdim)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 6L,  6L)
    current <- make_capped_volume_hypercube_box(8820L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(10289L, refdim)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  6L)
    current <- make_capped_volume_hypercube_box(10290L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(12004L, refdim)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 7L,  7L)
    current <- make_capped_volume_hypercube_box(12005L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(13719L, refdim)
    checkIdentical(target, current)

    target <- c( 7L,  7L, 5L, 8L,  7L)
    current <- make_capped_volume_hypercube_box(13720L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(15679L, refdim)
    checkIdentical(target, current)

    target <- c( 8L,  7L, 5L, 8L,  7L)
    current <- make_capped_volume_hypercube_box(15680L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(17919L, refdim)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  7L)
    current <- make_capped_volume_hypercube_box(17920L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(20479L, refdim)
    checkIdentical(target, current)

    target <- c( 8L,  8L, 5L, 8L,  8L)
    current <- make_capped_volume_hypercube_box(20480L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(23039L, refdim)
    checkIdentical(target, current)

    target <- c( 9L,  8L, 5L, 8L,  8L)
    current <- make_capped_volume_hypercube_box(23040L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(25919L, refdim)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  8L)
    current <- make_capped_volume_hypercube_box(25920L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(29159L, refdim)
    checkIdentical(target, current)

    target <- c( 9L,  9L, 5L, 8L,  9L)
    current <- make_capped_volume_hypercube_box(29160L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(32399L, refdim)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L,  9L)
    current <- make_capped_volume_hypercube_box(32400L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(35999L, refdim)
    checkIdentical(target, current)

    target <- c( 9L, 10L, 5L, 8L, 10L)
    current <- make_capped_volume_hypercube_box(36000L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(39999L, refdim)
    checkIdentical(target, current)

    target <- c( 10L, 10L, 5L, 8L, 10L)
    current <- make_capped_volume_hypercube_box(40000L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(43999L, refdim)
    checkIdentical(target, current)

    target <- c( 11L, 10L, 5L, 8L, 10L)
    current <- make_capped_volume_hypercube_box(44000L, refdim)
    checkIdentical(target, current)
    current <- make_capped_volume_hypercube_box(47999L, refdim)
    checkIdentical(target, current)

    refvol <- as.integer(prod(refdim))

    target <- c( 14L, 10L, 5L, 8L, 10L)
    current <- make_capped_volume_hypercube_box(refvol-1L, refdim)
    checkIdentical(target, current)

    current <- make_capped_volume_hypercube_box(refvol, refdim)
    checkIdentical(refdim, current)
}

