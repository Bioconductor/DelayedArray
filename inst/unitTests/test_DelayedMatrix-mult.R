#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

### Temporary workaround to make read_block(..., as.sparse=TRUE) work on
### a SVT_SparseArray object. Won't be needed anymore once read_block()
### gets switch from OLD_extract_sparse_array() to extract_sparse_array().
### TODO: Make sure to get rid of this when switching read_block()
### from OLD_extract_sparse_array() to extract_sparse_array().
setAs("COO_SparseArray", "SparseArraySeed",
    function(from)
        SparseArraySeed(dim(from), nzindex=nzcoo(from),
                                   nzdata=nzdata(from),
                                   dimnames=dimnames(from))
)
setMethod(OLD_extract_sparse_array, "SVT_SparseArray",
    function(x, index) {
        svt <- extract_sparse_array(x, index)
        as(as(svt, "COO_SparseArray"), "SparseArraySeed")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test BLOCK_mult_Lgrid() and BLOCK_mult_Rgrid()
###

test_BLOCK_mult_Lgrid_Rgrid <- function()
{
    ## BLOCK_mult_Lgrid() and BLOCK_mult_Rgrid() are the workhorses behind
    ## block matrix multiplication between a matrix-like object and an ordinary
    ## matrix (or other supported matrix-like object e.g. sparseMatrix or
    ## SparseMatrix).
    BLOCK_mult_Lgrid <- DelayedArray:::BLOCK_mult_Lgrid
    BLOCK_mult_Rgrid <- DelayedArray:::BLOCK_mult_Rgrid

    do_checks <- function(expected, x, y,
                          block_len, as.sparse, BPPARAM=NULL,
                          op=c("mult", "crossprod", "tcrossprod"),
                          BACKEND=NULL)
    {
        op <- match.arg(op)
        grid1 <- defaultAutoGrid(x, block.length=block_len)
        current <- BLOCK_mult_Lgrid(x, y, Lgrid=grid1, as.sparse=as.sparse,
                                          BPPARAM=BPPARAM, op=op,
                                          BACKEND=BACKEND)
        checkEquals(as.matrix(current), expected)
        grid2 <- defaultAutoGrid(y, block.length=block_len)
        current <- BLOCK_mult_Rgrid(x, y, Rgrid=grid2, as.sparse=as.sparse,
                                          BPPARAM=BPPARAM, op=op,
                                          BACKEND=BACKEND)
        checkEquals(as.matrix(current), expected)
    }

    ## Serial evaluation of:
    ##   <matrix> %*% <matrix>
    ##   <matrix> %*% <SVT_SparseMatrix>
    ##   <SVT_SparseMatrix> %*% <matrix>
    ##   <SVT_SparseMatrix> %*% <SVT_SparseMatrix>

    library(HDF5Array)
    m1 <- matrix(1:12, ncol=3, dimnames=list(letters[1:4], NULL))
    m2 <- matrix(101:115, nrow=3, dimnames=list(NULL, LETTERS[1:5]))
    m12 <- m1 %*% m2
    for (block_len in c(1:4, 6L, length(m1), length(m2), 1000L)) {
      for (x in list(m1, as(m1, "SVT_SparseMatrix"))) {
        for (y in list(m2, as(m2, "SVT_SparseMatrix"))) {
          do_checks(m12, x, y, block_len, as.sparse=FALSE)
          do_checks(m12, x, y, block_len, as.sparse=FALSE, BACKEND="HDF5Array")
          do_checks(m12, x, y, block_len, as.sparse=TRUE)
          do_checks(m12, x, y, block_len, as.sparse=TRUE, BACKEND="HDF5Array")
          do_checks(m12, t(x), y, block_len, as.sparse=FALSE, op="crossprod")
          do_checks(m12, t(x), y, block_len, as.sparse=FALSE, op="crossprod",
                                  BACKEND="HDF5Array")
          do_checks(m12, t(x), y, block_len, as.sparse=TRUE, op="crossprod")
          do_checks(m12, t(x), y, block_len, as.sparse=TRUE, op="crossprod",
                                  BACKEND="HDF5Array")
          do_checks(m12, x, t(y), block_len, as.sparse=FALSE, op="tcrossprod")
          do_checks(m12, x, t(y), block_len, as.sparse=FALSE, op="tcrossprod",
                                  BACKEND="HDF5Array")
          do_checks(m12, x, t(y), block_len, as.sparse=TRUE, op="tcrossprod")
          do_checks(m12, x, t(y), block_len, as.sparse=TRUE, op="tcrossprod",
                                  BACKEND="HDF5Array")
        }
      }
    }

    ## Parallel evaluation of:
    ##   <matrix> %*% <matrix>

    snow2 <- BiocParallel::SnowParam(workers=2)
    do_checks(m12, m1, m2, block_len=1L,
              as.sparse=FALSE, BPPARAM=snow2)
    do_checks(m12, m1, m2, block_len=1L,
              as.sparse=FALSE, BPPARAM=snow2, BACKEND="HDF5Array")
    do_checks(m12, m1, m2, block_len=6L,
              as.sparse=FALSE, BPPARAM=snow2)
    do_checks(m12, m1, m2, block_len=6L,
              as.sparse=FALSE, BPPARAM=snow2, BACKEND="HDF5Array")

    ## Serial evaluation of:
    ##   <DelayedMatrix> %*% <matrix>
    ##   <DelayedMatrix> %*% <SVT_SparseMatrix>
    ##   <matrix> %*% <DelayedMatrix>
    ##   <SVT_SparseMatrix> %*% <DelayedMatrix>

    M1 <- writeHDF5Array(m1, chunkdim=c(2, 2))
    M3 <- cbind(log(as(m2, "HDF5Array")), t(M1))
    m3 <- as.matrix(M3)
    m13 <- m1 %*% m3
    for (block_len in c(1:4, 6L, length(m1), length(m3), 1000L)) {
      for (y in list(m3, as(m3, "SVT_SparseMatrix"))) {
        grid1 <- defaultAutoGrid(M1, block.length=block_len)
        current <- BLOCK_mult_Lgrid(M1, y, Lgrid=grid1, as.sparse=FALSE,
                                           BPPARAM=NULL, BACKEND=NULL)
        checkEquals(current, m13)
        current <- BLOCK_mult_Lgrid(M1, y, Lgrid=grid1, as.sparse=FALSE,
                                           BPPARAM=NULL, BACKEND="HDF5Array")
        checkTrue(is(current, "DelayedMatrix"))
        checkTrue(validObject(current, complete=TRUE))
        checkEquals(as.matrix(current), m13)
      }
      for (x in list(m1, as(m1, "SVT_SparseMatrix"))) {
        grid3 <- defaultAutoGrid(M3, block.length=block_len)
        current <- BLOCK_mult_Rgrid(x, M3, Rgrid=grid3, as.sparse=FALSE,
                                           BPPARAM=NULL, BACKEND=NULL)
        checkEquals(current, m13)
        current <- BLOCK_mult_Rgrid(x, M3, Rgrid=grid3, as.sparse=FALSE,
                                           BPPARAM=NULL, BACKEND="HDF5Array")
        checkTrue(is(current, "DelayedMatrix"))
        checkTrue(validObject(current, complete=TRUE))
        checkEquals(as.matrix(current), m13)
      }
    }

    ## Parallel evaluation of:
    ##   <DelayedMatrix> %*% <matrix>
    ##   <matrix> %*% <DelayedMatrix>

    grid1 <- defaultAutoGrid(M1, block.length=1L)
    current <- BLOCK_mult_Lgrid(M1, m3, Lgrid=grid1, as.sparse=FALSE,
                                        BPPARAM=snow2, BACKEND=NULL)
    checkEquals(current, m13)
    current <- BLOCK_mult_Lgrid(M1, m3, Lgrid=grid1, as.sparse=FALSE,
                                        BPPARAM=snow2, BACKEND="HDF5Array")
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkEquals(as.matrix(current), m13)

    grid3 <- defaultAutoGrid(M3, block.length=2L)
    current <- BLOCK_mult_Rgrid(m1, M3, Rgrid=grid3, as.sparse=FALSE,
                                        BPPARAM=snow2, BACKEND=NULL)
    checkEquals(current, m13)
    current <- BLOCK_mult_Rgrid(m1, M3, Rgrid=grid3, as.sparse=FALSE,
                                        BPPARAM=snow2, BACKEND="HDF5Array")
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkEquals(as.matrix(current), m13)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test %*%, crossprod(), and tcrossprod() methods
###

TEST_m0 <- matrix(runif(60), ncol=6)
TEST_m0[2, 4] <- NA
TEST_m0[5, 4] <- Inf
TEST_m0[6, 3] <- -Inf
TEST_blocklengths0 <- c(3L, 5L, 15L, 1000L)
TEST_blocksizes0 <- TEST_blocklengths0 * get_type_size(type(TEST_m0))

### <matrix> %*% <DelayedMatrix> and <DelayedMatrix> %*% <matrix>
test_DelayedMatrix_mult <- function()
{
    M <- DelayedArray(realize(TEST_m0))

    Lm <- rbind(rep(1L, 10), rep(c(1L, 0L), 5), rep(-100L, 10))
    Rm <- rbind(Lm + 7.05, 0.1 * Lm)[ , 1:5]

    # Throws errors correctly.
    checkException(Lm[ , -1] %*% M, msg="non-conformable")

    ## With a dense DelayedMatrix.
    on.exit(suppressMessages(setAutoBlockSize()))
    for (block_size in TEST_blocksizes0) {
        suppressMessages(setAutoBlockSize(block_size))

        P <- Lm %*% M
        checkEquals(Lm %*% TEST_m0, as.matrix(P))
        P <- M %*% Rm
        checkEquals(TEST_m0 %*% Rm, as.matrix(P))

        P <- crossprod(t(Lm), M)
        checkEquals(crossprod(t(Lm), TEST_m0), as.matrix(P))
        P <- tcrossprod(M, t(Rm))
        checkEquals(tcrossprod(TEST_m0, t(Rm)), as.matrix(P))
    }

    ## With a sparse DelayedMatrix.
    s1 <- Matrix::rsparsematrix(nrow(TEST_m0), ncol(TEST_m0), density=0.2)
    s2 <- as(s1, "SVT_SparseMatrix")
    s0 <- as.matrix(s1)
    for (S in list(DelayedArray(s1), DelayedArray(s2))) {
      for (block_size in TEST_blocksizes0) {
        suppressMessages(setAutoBlockSize(block_size))

        P <- Lm %*% S
        checkEquals(unname(as.matrix(Lm %*% s0)), as.matrix(P))
        P <- S %*% Rm
        checkEquals(unname(as.matrix(s0 %*% Rm)), as.matrix(P))

        P <- crossprod(t(Lm), S)
        checkEquals(unname(as.matrix(crossprod(t(Lm), s0))), as.matrix(P))
        P <- tcrossprod(S, t(Rm))
        checkEquals(unname(as.matrix(tcrossprod(s0, t(Rm)))), as.matrix(P))
      }
    }

    ## Parallel evaluation.
    setAutoBPPARAM(BiocParallel::SnowParam(workers=2))
    on.exit(setAutoBPPARAM(), add=TRUE)
    suppressMessages(setAutoBlockSize(20))
    P <- Lm %*% M
    checkEquals(Lm %*% TEST_m0, as.matrix(P))
    P <- M %*% Rm
    checkEquals(TEST_m0 %*% Rm, as.matrix(P))
}

### <DelayedMatrix> %*% <DelayedMatrix>
### Based on DelayedArray:::.super_BLOCK_mult().
test_DelayedMatrix_mult_DelayedMatrix <- function()
{
    y1 <- matrix(runif(100), ncol=5)
    Y1 <- DelayedArray(realize(y1))
    y2 <- matrix(runif(100), nrow=5)
    Y2 <- DelayedArray(realize(y2))
    y12 <- y1 %*% y2

    for (block_size in TEST_blocksizes0) {
        suppressMessages(setAutoBlockSize(block_size))

        out <- Y1 %*% Y2
        if (!is.null(getAutoRealizationBackend()))
            out <- as.matrix(out)
        checkEquals(y12, out)

        out <- crossprod(t(Y1), Y2)
        if (!is.null(getAutoRealizationBackend()))
            out <- as.matrix(out)
        checkEquals(y12, out)

        out <- tcrossprod(Y1, t(Y2))
        if (!is.null(getAutoRealizationBackend()))
            out <- as.matrix(out)
        checkEquals(y12, out)
    }
}

### Based on DelayedArray:::.super_BLOCK_self()
test_DelayedMatrix_crossprod_self <- function()
{
    M <- DelayedArray(realize(TEST_m0))

    on.exit(DelayedArray:::setAutoMultParallelAgnostic())
    on.exit(suppressMessages(setAutoBlockSize()), add=TRUE)
    on.exit(setAutoBPPARAM(), add=TRUE)

    ## Check self-product in non-core-agnostic and core-agnostic ways.
    for (agnostic in c(FALSE, TRUE)) {
        ## Serial evaluation.
        setAutoBPPARAM()
        DelayedArray:::setAutoMultParallelAgnostic(agnostic)
        for (block_size in TEST_blocksizes0) {
            suppressMessages(setAutoBlockSize(block_size))
            checkEquals(as.matrix(crossprod(M)), crossprod(TEST_m0))
            checkEquals(as.matrix(tcrossprod(M)), tcrossprod(TEST_m0))
        }

        ## Parallel evaluation.
        setAutoBPPARAM(BiocParallel::SnowParam(workers=2))
        suppressMessages(setAutoBlockSize(20))
        checkEquals(as.matrix(crossprod(M)), crossprod(TEST_m0))
        checkEquals(as.matrix(tcrossprod(M)), tcrossprod(TEST_m0))
    }
}

