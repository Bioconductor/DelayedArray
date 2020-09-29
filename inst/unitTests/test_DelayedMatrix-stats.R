#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

a1 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))  # integer array
a2 <- a1 + runif(150) - 0.5                              # numeric array
m2 <- matrix(runif(60), ncol=6)                          # numeric matrix
m0 <- rsparsematrix(nrow=10, ncol=6, density=0.25)       # dgCMatrix object

block_sizes1 <- c(12L, 20L, 50L, 15000L)
block_sizes2 <- 2L * block_sizes1

test_DelayedMatrix_row_col_summarization <- function()
{
    test_row_col_summary <- function(FUN, m, M, block_sizes) {
        on.exit(suppressMessages(setAutoBlockSize()))
        FUN <- match.fun(FUN)
 
        target1 <- FUN(m)
        target2 <- FUN(m, na.rm=TRUE)
        target3 <- FUN(t(m))
        target4 <- FUN(t(m), na.rm=TRUE)
        for (block_size in block_sizes) {
            suppressMessages(setAutoBlockSize(block_size))
            current <- FUN(M)
            checkEquals(target1, current)
            checkIdentical(typeof(target1), typeof(current))
            current <- FUN(M, na.rm=TRUE)
            checkEquals(target2, current)
            checkIdentical(typeof(target2), typeof(current))
            current <- FUN(t(M))
            checkEquals(target3, current)
            checkIdentical(typeof(target3), typeof(current))
            current <- FUN(t(M), na.rm=TRUE)
            checkEquals(target4, current)
            checkIdentical(typeof(target4), typeof(current))
        }
    }

    FUNS <- c("rowSums", "colSums", "rowMeans", "colMeans",
              "rowMins", "colMins", "rowMaxs", "colMaxs",
              "rowRanges", "colRanges")

    ## on an integer matrix
    m <- a1[ , , 1]
    A1 <- realize(a1)
    M <- A1[ , , 1]
    for (FUN in FUNS) {
        test_row_col_summary(FUN, m, M, block_sizes2)
        test_row_col_summary(FUN, m[ , 0], M[ , 0], block_sizes2)
    }

    ## on a numeric matrix
    m <- m2
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    M <- realize(m)
    for (FUN in FUNS)
        test_row_col_summary(FUN, m, M, block_sizes2)
    library(genefilter)
    ## Note that the matrixStats package also defines a rowVars() function.
    test_row_col_summary(genefilter::rowVars, m, M, block_sizes2)

    ## on a sparse matrix
    m <- m0
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    M <- DelayedArray(m)
    for (FUN in FUNS)
        test_row_col_summary(FUN, as.matrix(m), M, block_sizes2)
    test_row_col_summary(genefilter::rowVars, m, M, block_sizes2)
}

