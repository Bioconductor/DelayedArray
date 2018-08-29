.TEST_MATRIX1 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6,
                   dimnames=list(NULL, LETTERS[1:6]))

.TEST_SAS2 <- SparseArraySeed(dim=c(40, 100, 1),
                              aind=arrayInd(1:50, c(10, 5, 1)),
                              nzdata=1:50)
.TEST_ARRAY2 <- sparse2dense(.TEST_SAS2)

test_DelayedSubset <- function()
{
    new_DelayedSubset <- DelayedArray:::new_DelayedSubset

    ## Ordinary array seed - no-op

    x1 <- new_DelayedSubset(.TEST_MATRIX1)
    a1 <- .TEST_MATRIX1

    checkTrue(is(x1, "DelayedSubset"))
    checkTrue(validObject(x1))
    checkIdentical(.TEST_MATRIX1, seed(x1))
    checkTrue(isNoOp(x1))
    checkIdentical(dim(a1), dim(x1))
    checkIdentical(dimnames(a1), dimnames(x1))
    checkIdentical(a1, as.array(x1))
    checkIdentical(a1[5:2, c(1:2, 5, 2:1)],
                   extract_array(x1, list(5:2, c(1:2, 5, 2:1))))
    checkIdentical(a1[5:2, 2L, drop=FALSE],
                   extract_array(x1, list(5:2, 2L)))
    checkIdentical(a1[5:2, integer(0)],
                   extract_array(x1, list(5:2, integer(0))))

    checkIdentical(FALSE, isSparse(x1))

    ## Ordinary array seed

    x2 <- new_DelayedSubset(.TEST_MATRIX1, list(5:2, c(1:2, 5, 2:1)))
    a2 <- .TEST_MATRIX1[5:2, c(1:2, 5L, 2:1)]

    checkTrue(is(x2, "DelayedSubset"))
    checkTrue(validObject(x2))
    checkIdentical(.TEST_MATRIX1, seed(x2))
    checkIdentical(FALSE, isNoOp(x2))
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))
    checkIdentical(a2[ , 4:2],
                   extract_array(x2, list(NULL, 4:2)))
    checkIdentical(a2[ , 4L, drop=FALSE],
                   extract_array(x2, list(NULL, 4L)))
    checkIdentical(a2[integer(0), ],
                   extract_array(x2, list(integer(0), NULL)))

    checkIdentical(FALSE, isSparse(x2))

    ## Sparse seed - no-op

    i1 <- c(3:8, 1L)
    x3 <- new_DelayedSubset(.TEST_SAS2)
    a3 <- .TEST_ARRAY2

    checkTrue(is(x3, "DelayedSubset"))
    checkTrue(validObject(x3))
    checkIdentical(.TEST_SAS2, seed(x3))
    checkTrue(isNoOp(x3))
    checkIdentical(dim(a3), dim(x3))
    checkIdentical(a3, as.array(x3))
    checkIdentical(a3[i1, 6:5, , drop=FALSE],
                   extract_array(x3, list(i1, 6:5, NULL)))
    checkIdentical(a3[i1, c(6:5, 6L), , drop=FALSE],
                   extract_array(x3, list(i1, c(6:5, 6L), NULL)))
    checkIdentical(a3[i1, c(6:5, 6L), integer(0)],
                   extract_array(x3, list(i1, c(6:5, 6L), integer(0))))

    checkTrue(isSparse(x3))
    sas3 <- .TEST_SAS2
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x3, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas3, current)
    current <- extract_sparse_array(x3, list(i1, 6:5, NULL))
    target <- extract_sparse_array(sas3, list(i1, 6:5, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x3, list(i1, 6:5, integer(0)))
    target <- extract_sparse_array(sas3, list(i1, 6:5, integer(0)))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)

    ## Sparse seed but structural sparsity NOT propagated because subscripts
    ## in 'index' contain duplicates

    x4 <- new_DelayedSubset(.TEST_SAS2, list(i1, c(6:5, 6L), NULL))
    a4 <- .TEST_ARRAY2[i1, c(6:5, 6L), , drop=FALSE]

    checkTrue(is(x4, "DelayedSubset"))
    checkTrue(validObject(x4))
    checkIdentical(.TEST_SAS2, seed(x4))
    checkIdentical(FALSE, isNoOp(x4))
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(a4, as.array(x4))

    checkIdentical(FALSE, isSparse(x4))

    ## Sparse seed -- structural sparsity propagated

    x5 <- new_DelayedSubset(.TEST_SAS2, list(i1, 6:5, NULL))
    a5 <- .TEST_ARRAY2[i1, 6:5, , drop=FALSE]

    checkTrue(is(x5, "DelayedSubset"))
    checkTrue(validObject(x5))
    checkIdentical(.TEST_SAS2, seed(x5))
    checkIdentical(FALSE, isNoOp(x5))
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(a5, as.array(x5))
    checkIdentical(a5[7:5, , , drop=FALSE],
                   extract_array(x5, list(7:5, NULL, NULL)))
    checkIdentical(a5[7:5, 2L, integer(0), drop=FALSE],
                   extract_array(x5, list(7:5, 2L, integer(0))))

    checkTrue(isSparse(x5))
    sas5 <- extract_sparse_array(.TEST_SAS2, list(i1, 6:5, NULL))
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x5, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas5, current)
    current <- extract_sparse_array(x5, list(c(6:7, 3:2), NULL, NULL))
    target <- extract_sparse_array(sas5, list(c(6:7, 3:2), NULL, NULL))
    checkIdentical(target, current)
    current <- extract_sparse_array(x5, list(c(6:7, 3:2), NULL, integer(0)))
    target <- extract_sparse_array(sas5, list(c(6:7, 3:2), NULL, integer(0)))
    checkIdentical(target, current)
    current <- extract_sparse_array(x5, list(c(6:7, 3:2), 2L, NULL))
    target <- extract_sparse_array(sas5, list(c(6:7, 3:2), 2L, NULL))
    checkIdentical(target, current)
}

