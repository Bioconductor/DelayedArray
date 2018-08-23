
test_SparseArray_getters <- function()
{
    aind1 <- aind2 <- rbind(c(2,4,3), c(2,1,3), c(5,4,3), c(5,3,3),
                            c(5,4,1), c(5,1,1), c(5,4,2), c(5,4,1))
    nzdata <- seq_len(nrow(aind1)) / 10
    sparse_array <- SparseArray(5:3, aind1, nzdata)

    checkIdentical(dim(sparse_array), 5:3)
    checkIdentical(length(sparse_array), 60L)
    storage.mode(aind2) <- "integer"
    checkIdentical(aind(sparse_array), aind2)
    checkIdentical(nzdata(sparse_array), nzdata)
}

test_dense2sparse_and_sparse2dense <- function()
{
    m <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
    sparse_array <- dense2sparse(m)
    checkTrue(is(sparse_array, "SparseArray"))
    checkIdentical(dim(sparse_array), dim(m))
    checkIdentical(sparse2dense(sparse_array), m)
}

