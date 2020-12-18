#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

### Toy integer 3D SparseArraySeed.
.make_toy_sas1 <- function()
{
    dim1 <- 5:3
    nzindex1 <- Lindex2Mindex(1:prod(dim1), dim1)
    set.seed(123)
    nzindex1 <- nzindex1[sample(nrow(nzindex1), 30), , drop=FALSE]
    nzdata1 <- 24:-5
    nzdata1[2:3] <- NA
    nzdata1[4:5] <- 0L
    SparseArraySeed(dim1, nzindex1, nzdata1)
}

### Toy integer 2D SparseArraySeed with no zeros or NAs.
.make_toy_sas1b <- function()
{
    dim1b <- 5:4
    nzindex1b <- Lindex2Mindex(1:prod(dim1b), dim1b)
    set.seed(123)
    nzindex1b <- nzindex1b[sample(nrow(nzindex1b)), , drop=FALSE]
    nzdata1b <- sample(10L, nrow(nzindex1b), replace=TRUE)
    SparseArraySeed(dim1b, nzindex1b, nzdata1b)
}

### Toy numeric 3D SparseArraySeed.
.make_toy_sas2 <- function()
{
    sas1 <- .make_toy_sas1()
    nzdata2 <- sas1@nzdata + 0.01
    nzdata2[2:3] <- NA
    nzdata2[4:5] <- 0
    nzdata2[6:8] <- c(NaN, Inf, -Inf)
    SparseArraySeed(dim(sas1), sas1@nzindex, nzdata2)
}

### Toy character 3D SparseArraySeed.
.make_toy_sas3 <- function()
{
    sas1 <- .make_toy_sas1()
    nzdata3 <- paste0(sas1@nzdata, "aXb")
    nzdata3[2:3] <- NA
    nzdata3[4:5] <- ""
    SparseArraySeed(dim(sas1), sas1@nzindex, nzdata3)
}

### Toy logical 3D SparseArraySeed.
.make_toy_sas4 <- function()
{
    sas1 <- .make_toy_sas1()
    set.seed(123)
    nzdata4 <- sample(c(TRUE, FALSE, NA), length(sas1@nzdata), replace=TRUE)
    SparseArraySeed(dim(sas1), sas1@nzindex, nzdata4)
}

test_abind_SparseArraySeed_objects <- function()
{
    abind_SparseArraySeed_objects <-
        DelayedArray:::abind_SparseArraySeed_objects
    simple_abind <- DelayedArray:::simple_abind

    test_on <- function(..., along) {
	objects <- list(...)
        current <- abind_SparseArraySeed_objects(objects, along)
        checkTrue(is(current, "SparseArraySeed"))
        arrays <- lapply(unname(objects), as.array)
        target <- do.call(simple_abind, c(arrays, list(along=along)))
        checkIdentical(target, as.array(current))
    }
    do_tests <- function(x1, x2, x3, x4, along) {
        test_on(x1, along=along)
        test_on(x2, along=along)
        test_on(x3, along=along)
        test_on(x4, along=along)
        test_on(x1, x2, along=along)
        test_on(x2, x1, along=along)
        test_on(x1, x3, along=along)
        test_on(x3, x1, along=along)
        test_on(x1, x4, along=along)
        test_on(x4, x1, along=along)
        test_on(x2, x3, along=along)
        test_on(x3, x2, along=along)
        test_on(x2, x4, along=along)
        test_on(x4, x2, along=along)
        test_on(x3, x4, along=along)
        test_on(x4, x3, along=along)
        test_on(x1, x2, x3, along=along)
        test_on(x1, x3, x2, along=along)
        test_on(x2, x1, x3, along=along)
        test_on(x2, x3, x1, along=along)
        test_on(x3, x1, x2, along=along)
        test_on(x3, x2, x1, along=along)
        test_on(x1, x2, x3, x4, along=along)
    }

    sas1 <- .make_toy_sas1()  # integer 3D SparseArraySeed
    sas2 <- .make_toy_sas2()  # numeric 3D SparseArraySeed
    dimnames(sas2) <- list(NULL, letters[1:4], letters[24:26])
    sas3 <- .make_toy_sas3()  # character 3D SparseArraySeed
    dimnames(sas3) <- list(NULL, LETTERS[1:4], NULL)
    sas4 <- .make_toy_sas4()  # logical 3D SparseArraySeed

    x1 <- extract_sparse_array(sas1, list(1L, NULL, NULL))
    x2 <- extract_sparse_array(sas1, list(1:3, NULL, NULL))
    x3 <- extract_sparse_array(sas1, list(integer(0), NULL, NULL))
    x4 <- extract_sparse_array(sas1, list(1:2, NULL, NULL))
    do_tests(x1, x2, x3, x4, along=1L)

    x1 <- extract_sparse_array(sas1, list(NULL, 1L, NULL))
    x2 <- extract_sparse_array(sas1, list(NULL, 1:3, NULL))
    x3 <- extract_sparse_array(sas1, list(NULL, integer(0), NULL))
    x4 <- extract_sparse_array(sas1, list(NULL, 1:2, NULL))
    do_tests(x1, x2, x3, x4, along=2L)

    x1 <- extract_sparse_array(sas1, list(NULL, NULL, 1L))
    x2 <- extract_sparse_array(sas1, list(NULL, NULL, 1:3))
    x3 <- extract_sparse_array(sas1, list(NULL, NULL, integer(0)))
    x4 <- extract_sparse_array(sas1, list(NULL, NULL, 1:2))
    do_tests(x1, x2, x3, x4, along=3L)
}

test_SparseArraySeed_unary_iso_ops <- function()
{
    do_tests <- function(.Generic, sas) {
        GENERIC <- match.fun(.Generic)
        current <- GENERIC(sas)
        checkTrue(is(current, "SparseArraySeed"))
        checkIdentical(dim(sas), dim(current))
        checkIdentical(GENERIC(as.array(sas)), as.array(current))
    }

    sas1 <- .make_toy_sas1()  # integer 3D SparseArraySeed
    sas2 <- .make_toy_sas2()  # numeric 3D SparseArraySeed
    sas4 <- .make_toy_sas4()  # logical 3D SparseArraySeed
    for (.Generic in c("is.na", "is.infinite", "is.nan")) {
        do_tests(.Generic, sas1)
        do_tests(.Generic, sas2)
        do_tests(.Generic, sas4)
    }

    sas3 <- .make_toy_sas3()  # character 3D SparseArraySeed
    for (.Generic in c("nchar", "tolower", "toupper")) {
        do_tests(.Generic, sas3)
    }
}

test_SparseArraySeed_anyNA <- function()
{
    do_test <- function(sas)
        checkIdentical(anyNA(as.array(sas)), anyNA(sas))

    sas1 <- .make_toy_sas1()    # integer 3D SparseArraySeed
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeros or NAs)
    sas2 <- .make_toy_sas2()    # numeric 3D SparseArraySeed
    sas3 <- .make_toy_sas3()    # character 3D SparseArraySeed
    sas4 <- .make_toy_sas4()    # logical 3D SparseArraySeed
    do_test(sas1)
    do_test(sas1b)
    do_test(sas2)
    do_test(sas3)
    do_test(sas4)
}

test_SparseArraySeed_which <- function()
{
    sas4 <- .make_toy_sas4()  # logical 3D SparseArraySeed
    checkIdentical(which(as.array(sas4)), which(sas4))
    checkIdentical(which(as.array(sas4), arr.ind=TRUE, useNames=FALSE),
                   which(sas4, arr.ind=TRUE))

    sas0 <- SparseArraySeed(5:4, nzindex=cbind(3:5, 4L), nzdata=FALSE)
    target1 <- integer(0)
    target2 <- matrix(integer(0), ncol=2)
    checkIdentical(target1, which(sas0))
    checkIdentical(target2, which(sas0, arr.ind=TRUE))
}

test_SparseArraySeed_Summary <- function()
{
    do_tests <- function(.Generic, sas) {
        GENERIC <- match.fun(.Generic)
        checkIdentical(GENERIC(as.array(sas)),
                       GENERIC(sas))
        checkIdentical(GENERIC(as.array(sas), na.rm=TRUE),
                       GENERIC(sas, na.rm=TRUE))
    }

    sas1 <- .make_toy_sas1()    # integer 3D SparseArraySeed
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeros or NAs)
    sas2 <- .make_toy_sas2()    # numeric 3D SparseArraySeed
    sas4 <- .make_toy_sas4()    # logical 3D SparseArraySeed
    for (.Generic in c("max", "min", "range", "sum", "prod")) {
        do_tests(.Generic, sas1)
        do_tests(.Generic, sas1b)
        do_tests(.Generic, sas2)
        do_tests(.Generic, sas4)
    }
    for (.Generic in c("any", "all")) {
        do_tests(.Generic, sas1)
        do_tests(.Generic, sas1b)
        suppressWarnings(do_tests(.Generic, sas2))
        do_tests(.Generic, sas4)
    }
}

test_SparseArraySeed_mean <- function()
{
    do_tests <- function(sas) {
        checkIdentical(mean(as.array(sas)), mean(sas))
        checkIdentical(mean(as.array(sas), na.rm=TRUE), mean(sas, na.rm=TRUE))
    }

    sas1 <- .make_toy_sas1()    # integer 3D SparseArraySeed
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeros or NAs)
    sas2 <- .make_toy_sas2()    # numeric 3D SparseArraySeed
    sas4 <- .make_toy_sas4()    # logical 3D SparseArraySeed
    do_tests(sas1)
    do_tests(sas1b)
    do_tests(sas2)
    do_tests(sas4)
}

