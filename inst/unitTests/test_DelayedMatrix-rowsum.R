#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

TEST_m1 <- matrix(runif(100), ncol=5, dimnames=list(NULL, LETTERS[1:5]))
TEST_m1[1, 2] <- NA
TEST_m1[2, 3] <- Inf
TEST_m1[3, 5] <- NaN
TEST_blocklengths1 <- c(1L, 6L, 10L, length(TEST_m1), 1000L * length(TEST_m1))
TEST_blocksizes1 <- TEST_blocklengths1 * get_type_size(type(TEST_m1))

TEST_group1 <- sample(6L, nrow(TEST_m1), replace=TRUE)

test_BLOCK_rowsum_colsum <- function()
{
    BLOCK_rowsum <- DelayedArray:::BLOCK_rowsum
    BLOCK_colsum <- DelayedArray:::BLOCK_colsum

    m1 <- TEST_m1
    m2 <- t(m1)
    group1 <- group2 <- TEST_group1

    library(HDF5Array)
    for (reorder in c(TRUE, FALSE)) {
        for (na.rm in c(FALSE, TRUE)) {
            rs1 <- rowsum(m1, group1, reorder=reorder, na.rm=na.rm)
            cs2 <- colsum(m2, group2, reorder=reorder, na.rm=na.rm)

            ## --- Serial evaluation ---

            for (block_len in TEST_blocklengths1) {
                grid1 <- defaultAutoGrid(m1, block.length=block_len)
                grid2 <- defaultAutoGrid(m2, block.length=block_len)

                ## In-memory realization.
                current <- BLOCK_rowsum(m1, group1,
                                        reorder=reorder, na.rm=na.rm,
                                        grid=grid1, as.sparse=NA,
                                        BPPARAM=NULL, BACKEND=NULL)
                checkEquals(current, rs1)
                current <- BLOCK_colsum(m2, group2,
                                        reorder=reorder, na.rm=na.rm,
                                        grid=grid2, as.sparse=NA,
                                        BPPARAM=NULL, BACKEND=NULL)
                checkEquals(current, cs2)

                ## On-disk realization (HDF5 file).
                current <- BLOCK_rowsum(m1, group1,
                                        reorder=reorder, na.rm=na.rm,
                                        grid=grid1, as.sparse=NA,
                                        BPPARAM=NULL, BACKEND="HDF5Array")
                checkTrue(is(current, "DelayedMatrix"))
                checkTrue(validObject(current, complete=TRUE))
                checkEquals(as.matrix(current), rs1)
                current <- BLOCK_colsum(m2, group2,
                                        reorder=reorder, na.rm=na.rm,
                                        grid=grid2, as.sparse=NA,
                                        BPPARAM=NULL, BACKEND="HDF5Array")
                checkTrue(is(current, "DelayedMatrix"))
                checkTrue(validObject(current, complete=TRUE))
                checkEquals(as.matrix(current), cs2)
            }

            ## --- Parallel evaluation ---

            snow2 <- BiocParallel::SnowParam(workers=2)
            grid1 <- defaultAutoGrid(m1, block.length=6L)
            grid2 <- defaultAutoGrid(m2, block.length=10L)

            ## In-memory realization.
            current <- BLOCK_rowsum(m1, group1, reorder=reorder, na.rm=na.rm,
                                    grid=grid1, as.sparse=NA,
                                    BPPARAM=snow2, BACKEND=NULL)
            checkEquals(current, rs1)
            current <- BLOCK_colsum(m2, group2, reorder=reorder, na.rm=na.rm,
                                    grid=grid2, as.sparse=NA,
                                    BPPARAM=snow2, BACKEND=NULL)
            checkEquals(current, cs2)

            ## On-disk realization (HDF5 file).
            current <- BLOCK_rowsum(m1, group1, reorder=reorder, na.rm=na.rm,
                                    grid=grid1, as.sparse=NA,
                                    BPPARAM=snow2, BACKEND="HDF5Array")
            checkTrue(is(current, "DelayedMatrix"))
            checkTrue(validObject(current, complete=TRUE))
            checkEquals(as.matrix(current), rs1)
            current <- BLOCK_colsum(m2, group2, reorder=reorder, na.rm=na.rm,
                                    grid=grid2, as.sparse=NA,
                                    BPPARAM=snow2, BACKEND="HDF5Array")
            checkTrue(is(current, "DelayedMatrix"))
            checkTrue(validObject(current, complete=TRUE))
            checkEquals(as.matrix(current), cs2)
        }
    }
}

test_rowsum_colsum_DelayedMatrix <- function()
{
    m1 <- TEST_m1
    m2 <- t(m1)
    group1 <- group2 <- TEST_group1

    M1 <- DelayedArray(realize(m1))
    M2 <- DelayedArray(realize(m2))

    on.exit(suppressMessages(setAutoBlockSize()))
    for (reorder in c(TRUE, FALSE)) {
        for (na.rm in c(FALSE, TRUE)) {
            rs1 <- rowsum(m1, group1, reorder=reorder, na.rm=na.rm)
            cs2 <- colsum(m2, group2, reorder=reorder, na.rm=na.rm)

            ## --- Serial evaluation ---

            for (block_size in TEST_blocksizes1) {
                suppressMessages(setAutoBlockSize(block_size))
                current <- rowsum(M1, group1, reorder=reorder, na.rm=na.rm)
                checkEquals(as.matrix(current), rs1)
                current <- colsum(M2, group2, reorder=reorder, na.rm=na.rm)
                checkEquals(as.matrix(current), cs2)
            }

            ## --- Parallel evaluation ---

            setAutoBPPARAM(BiocParallel::SnowParam(workers=2))
            on.exit(setAutoBPPARAM(), add=TRUE)
            suppressMessages(setAutoBlockSize(20))
            current <- rowsum(M1, group1, reorder=reorder, na.rm=na.rm)
            checkEquals(as.matrix(current), rs1)
            current <- colsum(M2, group2, reorder=reorder, na.rm=na.rm)
            checkEquals(as.matrix(current), cs2)
            setAutoBPPARAM()
       }
   }
}

