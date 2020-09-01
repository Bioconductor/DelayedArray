
test_dgCMatrix_stats <- function()
{
    colMins_dgCMatrix <- DelayedArray:::colMins_dgCMatrix
    colMaxs_dgCMatrix <- DelayedArray:::colMaxs_dgCMatrix
    colRanges_dgCMatrix <- DelayedArray:::colRanges_dgCMatrix
    colVars_dgCMatrix <- DelayedArray:::colVars_dgCMatrix

    do_other_tests <- function(m, m0) {
        ## colMins_dgCMatrix()
        target <- apply(m, 2, min)
        checkIdentical(target, colMins(m))
        checkIdentical(target, colMins_dgCMatrix(m0))
        target <- suppressWarnings(apply(m, 2, min, na.rm=TRUE))
        checkIdentical(target, colMins(m, na.rm=TRUE))
        checkIdentical(target, colMins_dgCMatrix(m0, na.rm=TRUE))

        ## colMaxs_dgCMatrix()
        target <- apply(m, 2, max)
        checkIdentical(target, colMaxs(m))
        checkIdentical(target, colMaxs_dgCMatrix(m0))
        target <- suppressWarnings(apply(m, 2, max, na.rm=TRUE))
        checkIdentical(target, colMaxs(m, na.rm=TRUE))
        checkIdentical(target, colMaxs_dgCMatrix(m0, na.rm=TRUE))

        ## colRanges_dgCMatrix()
        target <- t(apply(m, 2, range))
        checkIdentical(target, colRanges(m))
        checkIdentical(target, colRanges_dgCMatrix(m0))
        target <- suppressWarnings(t(apply(m, 2, range, na.rm=TRUE)))
        checkIdentical(target, colRanges(m, na.rm=TRUE))
        checkIdentical(target, colRanges_dgCMatrix(m0, na.rm=TRUE))

        ## colVars_dgCMatrix()
        target <- apply(m, 2, var)
        checkEquals(target, colVars(m))
        checkEquals(target, colVars_dgCMatrix(m0))
        target <- apply(m, 2, var, na.rm=TRUE)
        checkEquals(target, colVars(m, na.rm=TRUE))
        checkEquals(target, colVars_dgCMatrix(m0, na.rm=TRUE))
    }

    set.seed(123)
    m0 <- rsparsematrix(22, 10, density=0.25)
    m0[3L, 1L] <- NA
    m0[3L, 2L] <- NaN
    m0[4L, 4L] <- NA
    m0[8L, 4L] <- NaN
    m0[4L, 8L] <- NaN
    m0[22L, 8L] <- NA
    m0[ , 9L] <- NA
    m0[ , 10L] <- NaN
    m <- as.matrix(m0)

    ## rowsum()
    rgroup <- sample(9, nrow(m0), replace=TRUE)  # define groups of rows
    checkIdentical(rowsum(m, rgroup), rowsum(m0, rgroup))
    checkIdentical(rowsum(m, rgroup, reorder=FALSE),
                   rowsum(m0, rgroup, reorder=FALSE))
    checkIdentical(rowsum(m, rgroup, na.rm=TRUE),
                   rowsum(m0, rgroup, na.rm=TRUE))

    ## colsum()
    cgroup <- sample(5, ncol(m0), replace=TRUE)  # define groups of cols
    checkEquals(colsum(m, cgroup), colsum(m0, cgroup))
    checkEquals(colsum(m, cgroup, reorder=FALSE),
                colsum(m0, cgroup, reorder=FALSE))
    checkEquals(colsum(m, cgroup, na.rm=TRUE),
                colsum(m0, cgroup, na.rm=TRUE))

    do_other_tests(m, m0)
}

