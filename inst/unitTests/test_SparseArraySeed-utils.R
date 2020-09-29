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

### Toy integer 2D SparseArraySeed with no zeroes or NAs.
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
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeroes or NAs)
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
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeroes or NAs)
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
    sas1b <- .make_toy_sas1b()  # integer 2D SparseArraySeed (no zeroes or NAs)
    sas2 <- .make_toy_sas2()    # numeric 3D SparseArraySeed
    sas4 <- .make_toy_sas4()    # logical 3D SparseArraySeed
    do_tests(sas1)
    do_tests(sas1b)
    do_tests(sas2)
    do_tests(sas4)
}

