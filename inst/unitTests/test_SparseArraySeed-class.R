
test_SparseArraySeed_getters <- function()
{
    nzindex1 <- nzindex2 <- rbind(c(2,4,3), c(2,1,3), c(5,4,3), c(5,3,3),
                                  c(5,4,1), c(5,1,1), c(5,4,2), c(5,4,1))
    nzdata <- seq_len(nrow(nzindex1)) / 10
    sas <- SparseArraySeed(5:3, nzindex1, nzdata)

    checkIdentical(dim(sas), 5:3)
    checkIdentical(length(sas), 60L)
    storage.mode(nzindex2) <- "integer"
    checkIdentical(nzindex(sas), nzindex2)
    checkIdentical(nzdata(sas), nzdata)
}

test_dense2sparse_and_sparse2dense <- function()
{
    m <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
    sas <- dense2sparse(m)
    checkTrue(is(sas, "SparseArraySeed"))
    checkIdentical(dim(sas), dim(m))
    checkIdentical(sparse2dense(sas), m)
}

