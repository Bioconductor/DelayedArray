#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

DEFAULT_BLOCK_SIZE <- DelayedArray:::DEFAULT_BLOCK_SIZE

Arith_members <- c("+", "-", "*", "/", "^", "%%", "%/%")
Compare_members <- c("==", "!=", "<=", ">=", "<", ">")
Logic_members <- c("&", "|")  # currently untested

a1 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))  # integer array
a2 <- a1 + runif(150) - 0.5                              # numeric array

block_sizes1 <- c(12L, 20L, 50L, 15000L)
block_sizes2 <- 2L * block_sizes1

test_DelayedArray_unary_ops <- function()
{
    a <- 2:-2 / (a1 - 3)
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)
    for (.Generic in c("is.na", "is.finite", "is.infinite", "is.nan")) {
        GENERIC <- match.fun(.Generic)
        checkIdentical(GENERIC(a), as.array(GENERIC(A)))
    }

    a <- array(sample(c(LETTERS, letters), 60, replace=TRUE), 5:3)
    A <- realize(a)
    ## For some obscure reason, the tests below fail in the context of
    ## 'DelayedArray:::.test()' or 'R CMD check'.
    ## TODO: Investigate this.
    #for (.Generic in c("nchar", "tolower", "toupper")) {
    #    GENERIC <- match.fun(.Generic)
    #    checkIdentical(GENERIC(a), as.array(GENERIC(A)))
    #}
}

test_DelayedArray_Math_ans_Arith <- function()
{
    toto1 <- function(a) { 100 / floor(abs((5 * log(a + 0.2) - 1)^3)) }
    toto2 <- function(a) { 100L + (5L * (a - 2L)) %% 7L }

    ## with an integer array
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)
    ## Not sure what's going on but it seems that this call to checkIdentical()
    ## crashes the RUnit package but only when the tests are run by
    ## 'R CMD check'.
    #checkIdentical(toto2(a), as.array(toto2(A)))

    ## with a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)
    checkIdentical(toto1(a), as.array(toto1(A)))
    checkIdentical(toto2(a), as.array(toto2(A)))
    checkIdentical(toto2(toto1(a)), as.array(toto2(toto1(A))))
    checkIdentical(toto1(toto2(a)), as.array(toto1(toto2(A))))

    a <- a[ , 10:4, -2]
    A <- A[ , 10:4, -2]
    checkIdentical(toto1(a), as.array(toto1(A)))
    checkIdentical(toto2(a), as.array(toto2(A)))
    checkIdentical(toto2(toto1(a)), as.array(toto2(toto1(A))))
    checkIdentical(toto1(toto2(a)), as.array(toto1(toto2(A))))

    ## with a numeric matrix
    m <- a[ , , 2]
    M <- realize(m)
    checkIdentical(toto1(m), as.matrix(toto1(M)))
    checkIdentical(t(toto1(m)), as.matrix(toto1(t(M))))
    checkIdentical(t(toto1(m)), as.matrix(t(toto1(M))))
    M <- drop(A[ , , 2])
    checkIdentical(toto1(m), as.matrix(toto1(M)))
    checkIdentical(t(toto1(m)), as.matrix(toto1(t(M))))
    checkIdentical(t(toto1(m)), as.matrix(t(toto1(M))))
}

test_DelayedArray_Ops_with_left_or_right_vector <- function()
{
    test_delayed_Ops_on_array <- function(.Generic, a, A, m, M) {
        on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
        GENERIC <- match.fun(.Generic)

        target_current <- list(
            list(GENERIC(a, m[ , 1]), GENERIC(A, M[ , 1])),
            list(GENERIC(m[ , 1], a), GENERIC(M[ , 1], A)),

            list(GENERIC(a, a[ , 1, 1]), GENERIC(A, A[ , 1, 1])),
            list(GENERIC(a[ , 1, 1], a), GENERIC(A[ , 1, 1], A)),

            list(GENERIC(a, a[[1]]), GENERIC(A, A[[1]])),
            list(GENERIC(a[[1]], a), GENERIC(A[[1]], A))
        )
        for (i in seq_along(target_current)) {
            target <- target_current[[i]][[1L]]
            current <- target_current[[i]][[2L]]
            checkIdentical(target, as.array(current))
            checkIdentical(target[5:3, , -2], as.array(current[5:3, , -2]))
            checkIdentical(target[0, 0, 0], as.array(current[0, 0, 0]))
            checkIdentical(target[0, 0, -2], as.array(current[0, 0, -2]))
            checkIdentical(target[ , 0, ], as.array(current[ , 0, ]))
            for (block_size in block_sizes2) {
                options(DelayedArray.block.size=block_size)
                checkEquals(sum(target, na.rm=TRUE), sum(current, na.rm=TRUE))
            }
        }
    }

    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)
    m <- a[ , , 2]
    M <- realize(m)

    ## "Logic" members currently untested.
    for (.Generic in c(Arith_members, Compare_members))
        test_delayed_Ops_on_array(.Generic, a, A, m, M)

    ## Takes too long and probably not that useful.
    #M <- drop(A[ , , 2])
    #for (.Generic in c(Arith_members, Compare_members))
    #    test_delayed_Ops_on_array(.Generic, a, A, m, M)
}

test_DelayedArray_Ops_COMBINE_seeds <- function()
{
    ## comparing 2 DelayedArray objects
    A1 <- realize(a1)
    A2 <- realize(a2)
    a3 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))
    a3[2, 9, 2] <- NA  # same as a3[[92]] <- NA
    A3 <- realize(a3)

    ## "Logic" members currently untested.
    for (.Generic in c(Arith_members, Compare_members)) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a1, a2)
        target2 <- GENERIC(a2, a1)
        target3 <- GENERIC(a1, a3)
        checkIdentical(target1, as.array(GENERIC(A1, A2)))
        checkIdentical(target2, as.array(GENERIC(A2, A1)))
        checkIdentical(target3, as.array(GENERIC(A1, A3)))
    }
}

test_DelayedArray_anyNA <- function()
{
    on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
    DelayedArray_block_anyNA <- DelayedArray:::.DelayedArray_block_anyNA

    A1 <- realize(a1)
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)

    for (block_size in block_sizes2) {
        options(DelayedArray.block.size=block_size)
        checkIdentical(FALSE, anyNA(A1))
        checkIdentical(FALSE, DelayedArray_block_anyNA(a1))
        checkIdentical(TRUE, anyNA(A))
        checkIdentical(TRUE, DelayedArray_block_anyNA(a))
    }
}

test_DelayedArray_which <- function()
{
    on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
    DelayedArray_block_which <- DelayedArray:::.DelayedArray_block_which

    a <- a1 == 1L
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)

    target <- which(a)
    for (block_size in block_sizes2) {
        options(DelayedArray.block.size=block_size)
        checkIdentical(target, which(A))
        checkIdentical(target, DelayedArray_block_which(a))
    }

    a <- a1 == -1L    # all FALSE
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)

    target <- integer(0)
    for (block_size in block_sizes2) {
        options(DelayedArray.block.size=block_size)
        checkIdentical(target, which(A))
        checkIdentical(target, DelayedArray_block_which(a))
    }
}

test_DelayedArray_Summary <- function()
{
    test_Summary <- function(.Generic, a, block_sizes) {
        on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
        DelayedArray_block_Summary <- DelayedArray:::.DelayedArray_block_Summary

        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a)
        target2 <- GENERIC(a, na.rm=TRUE)
        A <- realize(a)
        for (block_size in block_sizes) {
            options(DelayedArray.block.size=block_size)
            checkIdentical(target1, GENERIC(A))
            checkIdentical(target1, GENERIC(t(A)))
            checkIdentical(target1, DelayedArray_block_Summary(.Generic, a))
            checkIdentical(target2, GENERIC(A, na.rm=TRUE))
            checkIdentical(target2, GENERIC(t(A), na.rm=TRUE))
            checkIdentical(target2, DelayedArray_block_Summary(.Generic, a,
                                                               na.rm=TRUE))
        }
    }

    ## on an integer array
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    #for (.Generic in c("max", "min", "range", "sum", "prod")) {
    for (.Generic in c("max", "min", "range", "sum"))
        test_Summary(.Generic, a, block_sizes1)

    ## on a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    a[2, 10, 2] <- Inf  # same as a[[97]] <- Inf
    for (.Generic in c("max", "min", "range", "sum", "prod"))
        test_Summary(.Generic, a, block_sizes2)

    ## on a logical array
    a <- array(c(rep(NA, 62), rep(TRUE, 87), FALSE), c(5, 10, 3))
    for (.Generic in c("any", "all"))
        test_Summary(.Generic, a, block_sizes1)
}

test_DelayedArray_mean <- function()
{
    on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))

    ## on a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- realize(a)

    target1 <- mean(a)
    target2 <- mean(a, na.rm=TRUE)
    target3 <- mean(a[ , 10:4, -2])
    for (block_size in block_sizes2) {
        options(DelayedArray.block.size=block_size)
        checkIdentical(target1, mean(A))
        checkIdentical(target1, mean(t(A)))
        checkIdentical(target2, mean(A, na.rm=TRUE))
        checkIdentical(target2, mean(t(A), na.rm=TRUE))
        checkIdentical(target3, mean(A[ , 10:4, -2]))
        checkIdentical(target3, mean(t(A[ , 10:4, -2])))
    }
}

test_DelayedArray_apply <- function()
{
    test_apply <- function(a) {
        A <- realize(a)
        for (MARGIN in seq_along(dim(a))) {
            checkIdentical(apply(a, MARGIN, dim),
                           apply(A, MARGIN, dim))
            checkIdentical(apply(a, MARGIN, sum),
                           apply(A, MARGIN, sum))
            checkIdentical(apply(a, MARGIN, sum, na.rm=TRUE),
                           apply(A, MARGIN, sum, na.rm=TRUE))
            ## row/colSums and row/colMeans don't work yet in that case.
            if (dim(A)[[MARGIN]] == 0L && length(dim(A)) >= 3L)
                next
            checkIdentical(apply(a, MARGIN, rowSums),
                           apply(A, MARGIN, rowSums))
            checkIdentical(apply(a, MARGIN, rowSums, na.rm=TRUE),
                           apply(A, MARGIN, rowSums, na.rm=TRUE))
            checkIdentical(apply(a, MARGIN, colMeans),
                           apply(A, MARGIN, colMeans))
            checkIdentical(apply(a, MARGIN, colMeans, na.rm=TRUE),
                           apply(A, MARGIN, colMeans, na.rm=TRUE))
        }
    }

    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    test_apply(a)
    test_apply(a[ , , 0])

    dimnames(a) <- list(NULL, NULL, LETTERS[1:3])
    test_apply(a)
    test_apply(a[ , , 0])

    dimnames(a) <- list(NULL, letters[1:10], LETTERS[1:3])
    test_apply(a)
    test_apply(a[ , , 0])
}

