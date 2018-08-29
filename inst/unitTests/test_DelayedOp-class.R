new_DelayedSubset <- DelayedArray:::new_DelayedSubset
new_DelayedAperm <- DelayedArray:::new_DelayedAperm


.TEST_MATRIX1 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6,
                        dimnames=list(NULL, LETTERS[1:6]))

.TEST_SAS2 <- SparseArraySeed(dim=c(40, 100, 1),
                              aind=arrayInd(1:50, c(10, 5, 1)),
                              nzdata=1:50)
.TEST_ARRAY2 <- sparse2dense(.TEST_SAS2)

.check_DelayedUnaryOp_with_MATRIX1_seed_is_no_op <- function(x)
{
    a <- .TEST_MATRIX1
    checkIdentical(dim(a), dim(x))
    checkIdentical(dimnames(a), dimnames(x))
    checkIdentical(a, as.array(x))
    checkIdentical(unname(a[5:2, c(1:2, 5, 2:1)]),
                   unname(extract_array(x, list(5:2, c(1:2, 5, 2:1)))))
    checkIdentical(unname(a[5:2, 2L, drop=FALSE]),
                   unname(extract_array(x, list(5:2, 2L))))
    checkIdentical(unname(a[5:2, integer(0)]),
                   unname(extract_array(x, list(5:2, integer(0)))))
}

.check_DelayedUnaryOp_with_SAS2_seed_is_no_op <- function(x)
{
    a <- .TEST_ARRAY2
    checkIdentical(dim(a), dim(x))
    checkIdentical(dimnames(a), dimnames(x))
    checkIdentical(a, as.array(x))
    i1 <- c(3:8, 1L)
    checkIdentical(unname(a[i1, 6:5, , drop=FALSE]),
                   unname(extract_array(x, list(i1, 6:5, NULL))))
    checkIdentical(unname(a[i1, c(6:5, 6L), , drop=FALSE]),
                   unname(extract_array(x, list(i1, c(6:5, 6L), NULL))))
    checkIdentical(unname(a[i1, c(6:5, 6L), integer(0)]),
                   unname(extract_array(x, list(i1, c(6:5, 6L), integer(0)))))

    checkTrue(isSparse(x))
    sas <- .TEST_SAS2
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas, current)
    current <- extract_sparse_array(x, list(i1, 6:5, NULL))
    target <- extract_sparse_array(sas, list(i1, 6:5, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x, list(i1, 6:5, integer(0)))
    target <- extract_sparse_array(sas, list(i1, 6:5, integer(0)))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
}

test_DelayedSubset_constructor <- function(silent=FALSE)
{
    ## We also test seed() and isNoOp()

    x <- new_DelayedSubset()
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(NULL), x@index)
    checkTrue(isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1)
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(NULL, 5:4))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, 5:4), x@index)
    checkIdentical(FALSE, isNoOp(x))

    ## Test normalization of user-supplied 'Nindex' argument

    checkException(new_DelayedSubset(.TEST_MATRIX1, 5:4),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX1, list(5:4)),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX1, list(NULL, 7:4)),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX1, list(NULL, "zz")),
                   silent=silent)

    x <- new_DelayedSubset(.TEST_MATRIX1, list(0, c(5:4, 5.99, 6.99)))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(integer(0), c(5:4, 5:6)), x@index)
    checkIdentical(FALSE, isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(NULL, -3))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, c(1:2, 4:6)), x@index)
    checkIdentical(FALSE, isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(NULL, 1:6))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(-22, TRUE))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(NULL, c("E", "B", "C", "C")))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(NULL, c(5L, 2L, 3L, 3L)), x@index)
    checkIdentical(FALSE, isNoOp(x))

    x <- new_DelayedSubset(.TEST_MATRIX1, list(IRanges(3:2, 5), IRanges(1, 6)))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(c(3:5, 2:5), NULL), x@index)
    checkIdentical(FALSE, isNoOp(x))

    checkException(new_DelayedSubset(.TEST_MATRIX1, list(NULL, IRanges(0, 3))),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX1, list(NULL, IRanges(2, 7))),
                   silent=silent)
}

test_DelayedSubset_API <- function()
{
    ## 1. Ordinary array seed - no-op

    x1 <- new_DelayedSubset(.TEST_MATRIX1)
    .check_DelayedUnaryOp_with_MATRIX1_seed_is_no_op(x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed

    Nindex2 <- list(5:2, c(1:2, 5, 2:1))
    x2 <- new_DelayedSubset(.TEST_MATRIX1, Nindex2)

    a2 <- .TEST_MATRIX1[5:2, c(1:2, 5L, 2:1)]
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))
    checkIdentical(unname(a2[ , 4:2]),
                   unname(extract_array(x2, list(NULL, 4:2))))
    checkIdentical(unname(a2[ , 4L, drop=FALSE]),
                   unname(extract_array(x2, list(NULL, 4L))))
    checkIdentical(unname(a2[integer(0), ]),
                   unname(extract_array(x2, list(integer(0), NULL))))

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed - no-op

    x3 <- new_DelayedSubset(.TEST_SAS2)
    .check_DelayedUnaryOp_with_SAS2_seed_is_no_op(x3)

    ## 4. Sparse seed -- structural sparsity propagated

    i1 <- c(3:8, 1L)
    Nindex4 <- list(i1, 6:5, NULL)
    x4 <- new_DelayedSubset(.TEST_SAS2, Nindex4)

    a4 <- .TEST_ARRAY2[i1, 6:5, , drop=FALSE]
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(dimnames(a4), dimnames(x4))
    checkIdentical(a4, as.array(x4))
    checkIdentical(unname(a4[7:5, , , drop=FALSE]),
                   unname(extract_array(x4, list(7:5, NULL, NULL))))
    checkIdentical(unname(a4[7:5, 2L, integer(0), drop=FALSE]),
                   unname(extract_array(x4, list(7:5, 2L, integer(0)))))

    checkTrue(isSparse(x4))
    sas4 <- extract_sparse_array(.TEST_SAS2, Nindex4)
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x4, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas4, current)
    current <- extract_sparse_array(x4, list(c(6:7, 3:2), NULL, NULL))
    target <- extract_sparse_array(sas4, list(c(6:7, 3:2), NULL, NULL))
    checkIdentical(target, current)
    current <- extract_sparse_array(x4, list(c(6:7, 3:2), NULL, integer(0)))
    target <- extract_sparse_array(sas4, list(c(6:7, 3:2), NULL, integer(0)))
    checkIdentical(target, current)
    current <- extract_sparse_array(x4, list(c(6:7, 3:2), 2L, NULL))
    target <- extract_sparse_array(sas4, list(c(6:7, 3:2), 2L, NULL))
    checkIdentical(target, current)

    ## 5. Sparse seed but structural sparsity NOT propagated because
    ##    subscripts in 'Nindex' argument contain duplicates

    i1 <- c(3:8, 1L)
    Nindex5 <- list(i1, c(6:5, 6L), NULL)  # Nindex5[[2]] contains duplicates!
    x5 <- new_DelayedSubset(.TEST_SAS2, Nindex5)

    a5 <- .TEST_ARRAY2[i1, c(6:5, 6L), , drop=FALSE]
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))

    checkIdentical(FALSE, isSparse(x5))  # structural sparsity not propagated!
}

test_DelayedAperm_constructor <- function(silent=FALSE)
{
    ## We also test seed() and isNoOp()

    x <- new_DelayedAperm()
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(1L, x@perm)
    checkTrue(isNoOp(x))

    x <- new_DelayedAperm(.TEST_SAS2)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(1:3, x@perm)
    checkTrue(isNoOp(x))

    x <- new_DelayedAperm(.TEST_SAS2, 3:1)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(3:1, x@perm)
    checkIdentical(FALSE, isNoOp(x))

    checkException(new_DelayedAperm(.TEST_SAS2, "2"), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS2, integer(0)), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS2, NA_integer_), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS2, 9L), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS2, 1L), silent=silent)

    x <- new_DelayedAperm(.TEST_SAS2, 2:1)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(2:1, x@perm)
    checkIdentical(FALSE, isNoOp(x))
}

test_DelayedAperm_API <- function()
{
    ## 1. Ordinary array seed - no-op

    x1 <- new_DelayedAperm(.TEST_MATRIX1)
    .check_DelayedUnaryOp_with_MATRIX1_seed_is_no_op(x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed - transposition

    x2 <- new_DelayedAperm(.TEST_MATRIX1, 2:1)

    a2 <- t(.TEST_MATRIX1)
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))
    checkIdentical(unname(a2[4:2, ]),
                   unname(extract_array(x2, list(4:2, NULL))))
    checkIdentical(unname(a2[4L, , drop=FALSE]),
                   unname(extract_array(x2, list(4L, NULL))))
    checkIdentical(unname(a2[ , integer(0)]),
                   unname(extract_array(x2, list(NULL, integer(0)))))

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedAperm(.TEST_SAS2)
    .check_DelayedUnaryOp_with_SAS2_seed_is_no_op(x3)

    ## 4. Sparse seed -- transpose 1st and 3rd dims

    x4 <- new_DelayedAperm(.TEST_SAS2, 3:1)

    a4 <- aperm(.TEST_ARRAY2, 3:1)
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(dimnames(a4), dimnames(x4))
    checkIdentical(a4, as.array(x4))
    i1 <- c(3:8, 1L)
    checkIdentical(unname(a4[ , 6:5, i1, drop=FALSE]),
                   unname(extract_array(x4, list(NULL, 6:5, i1))))
    checkIdentical(unname(a4[ , c(6:5, 6L), i1, drop=FALSE]),
                   unname(extract_array(x4, list(NULL, c(6:5, 6L), i1))))
    checkIdentical(unname(a4[integer(0), c(6:5, 6L), i1]),
                   unname(extract_array(x4, list(integer(0), c(6:5, 6L), i1))))

    checkTrue(isSparse(x4))
    sas4 <- aperm(.TEST_SAS2, 3:1)
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x4, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas4, current)
    current <- extract_sparse_array(x4, list(NULL, 6:5, i1))
    target <- extract_sparse_array(sas4, list(NULL, 6:5, i1))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x4, list(integer(0), 6:5, i1))
    target <- extract_sparse_array(sas4, list(integer(0), 6:5, i1))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)

    ## 5. Sparse seed -- transpose 1st and 2nd dims and drop 3rd dim

    x5 <- new_DelayedAperm(.TEST_SAS2, 2:1)

    a5 <- t(drop(.TEST_ARRAY2))
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))
    i1 <- c(3:8, 1L)
    checkIdentical(unname(a5[6:5, i1]),
                   unname(extract_array(x5, list(6:5, i1))))
    checkIdentical(unname(a5[c(6:5, 6L), i1]),
                   unname(extract_array(x5, list(c(6:5, 6L), i1))))
    checkIdentical(unname(a5[c(6:5, 6L), integer(0)]),
                   unname(extract_array(x5, list(c(6:5, 6L), integer(0)))))

    checkTrue(isSparse(x5))
    sas5 <- aperm(.TEST_SAS2, 2:1)
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x5, list(NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas5, current)
    current <- extract_sparse_array(x5, list(6:5, i1))
    target <- extract_sparse_array(sas5, list(6:5, i1))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x5, list(6:5, integer(0)))
    target <- extract_sparse_array(sas5, list(6:5, integer(0)))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
}

