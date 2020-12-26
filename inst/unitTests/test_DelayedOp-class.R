new_DelayedSubset <- DelayedArray:::new_DelayedSubset
new_DelayedAperm <- DelayedArray:::new_DelayedAperm
new_DelayedUnaryIsoOpStack <- DelayedArray:::new_DelayedUnaryIsoOpStack
new_DelayedUnaryIsoOpWithArgs <- DelayedArray:::new_DelayedUnaryIsoOpWithArgs
new_DelayedSetDimnames <- DelayedArray:::new_DelayedSetDimnames
.INHERIT_FROM_SEED <- DelayedArray:::.INHERIT_FROM_SEED
new_DelayedNaryIsoOp <- DelayedArray:::new_DelayedNaryIsoOp

.TEST_DIM2 <- c(5L, 6L)

.TEST_MATRIX2a <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6,
                         dimnames=list(NULL, LETTERS[1:6]))

.TEST_SAS2b <- SparseArraySeed(
    dim=.TEST_DIM2,
    nzindex=Lindex2Mindex(c(2:3, 12:22, 30), .TEST_DIM2),
    nzdata=rep(TRUE, 14))

.TEST_ARRAY2b <- sparse2dense(.TEST_SAS2b)

.TEST_SAS2c <- SparseArraySeed(
    dim=.TEST_DIM2,
    nzindex=Lindex2Mindex(c(8:16, 21), .TEST_DIM2),
    nzdata=rep(TRUE, 10))

.TEST_ARRAY2c <- sparse2dense(.TEST_SAS2c)

.TEST_SAS3 <- SparseArraySeed(
    dim=c(40, 100, 1),
    nzindex=Lindex2Mindex(1:50, c(10, 5, 1)),
    nzdata=1:50)

.TEST_ARRAY3 <- sparse2dense(.TEST_SAS3)


.TEST_DIM4 <- c(6L, 10L, 3L, 2L)

.TEST_ARRAY4a <- array(seq_len(prod(.TEST_DIM4)), .TEST_DIM4,
                       dimnames=list(NULL,
                                     NULL,
                                     c("z1", "z2", "z3"),
                                     c("a1", "a2")))

set.seed(99L)
.TEST_ARRAY4b <- array(runif(prod(.TEST_DIM4)), .TEST_DIM4,
                       dimnames=list(NULL,
                                     LETTERS[1:10],
                                     NULL,
                                     c("b1", "b2")))

.TEST_SAS4c <- SparseArraySeed(
    dim=.TEST_DIM4,
    nzindex=Lindex2Mindex(51:130, .TEST_DIM4),
    nzdata=51:130)
.TEST_ARRAY4c <- sparse2dense(.TEST_SAS4c)

.TEST_SAS4d <- SparseArraySeed(
    dim=.TEST_DIM4,
    nzindex=Lindex2Mindex(11:110, .TEST_DIM4),
    nzdata=runif(100, max=150))
.TEST_ARRAY4d <- sparse2dense(.TEST_SAS4d)

.basic_checks_on_DelayedOp_with_DIM2 <- function(a, x)
{
    ## We use suppressWarnings() to suppress the warnings that some
    ## calls to extract_array() could generate on some particular
    ## DelayedOp objects. For example on a DelayedUnaryIsoOpStack object
    ## with the log() function in its OPS stack and negative array elements
    ## in its seed. These warnings are expected but annoying in the context
    ## of unit tests.
    checkIdentical(dim(a), dim(x))
    checkIdentical(dimnames(a), dimnames(x))
    checkIdentical(a, suppressWarnings(as.array(x)))

    i2 <- c(1:2, 6:2, 2:1)
    current <- suppressWarnings(extract_array(x, list(5:2, i2)))
    checkIdentical(unname(a[5:2, i2]), unname(current))
    current <- suppressWarnings(extract_array(x, list(5:2, NULL)))
    checkIdentical(unname(a[5:2, ]), unname(current))
    current <- suppressWarnings(extract_array(x, list(5:2, 2L)))
    checkIdentical(unname(a[5:2, 2L, drop=FALSE]), unname(current))
    current <- suppressWarnings(extract_array(x, list(5:2, integer(0))))
    checkIdentical(unname(a[5:2, integer(0)]), unname(current))
}

.basic_checks_on_DelayedOp_with_DIM3 <- function(a, x)
{
    checkIdentical(dim(a), dim(x))
    checkIdentical(dimnames(a), dimnames(x))
    checkIdentical(a, as.array(x))

    i3 <- c(3:8, 1L)
    current <- extract_array(x, list(i3, 6:5, NULL))
    checkIdentical(unname(a[i3, 6:5, , drop=FALSE]), unname(current))
    current <- extract_array(x, list(i3, c(6:5, 6L), NULL))
    checkIdentical(unname(a[i3, c(6:5, 6L), , drop=FALSE]), unname(current))
    current <- extract_array(x, list(i3, c(6:5, 6L), integer(0)))
    checkIdentical(unname(a[i3, c(6:5, 6L), integer(0)]), unname(current))
}

.check_extract_sparse_array_on_DelayedOp_with_DIM3 <- function(a, x)
{
    checkTrue(is_sparse(x))

    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(a, sparse2dense(current))

    i3 <- c(3:8, 1L)
    current <- extract_sparse_array(x, list(i3, 6:5, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    target <- extract_array(a, list(i3, 6:5, NULL))
    checkIdentical(target, sparse2dense(current))

    current <- extract_sparse_array(x, list(i3, 6:5, integer(0)))
    checkTrue(is(current, "SparseArraySeed"))
    target <- extract_array(a, list(i3, 6:5, integer(0)))
    checkIdentical(target, sparse2dense(current))
}

.basic_checks_on_DelayedOp_with_DIM4 <- function(a, x)
{
    checkIdentical(dim(a), dim(x))
    checkIdentical(dimnames(a), dimnames(x))
    checkIdentical(a, as.array(x))

    i1 <- c(3:6, 1L)
    current <- extract_array(x, list(i1, 6:5, NULL, 2:1))
    checkIdentical(unname(a[i1, 6:5, , 2:1]),
                   unname(current))
    current <- extract_array(x, list(i1, c(6:5, 6L), NULL, 2:1))
    checkIdentical(unname(a[i1, c(6:5, 6L), , 2:1]),
                   unname(current))
    current <- extract_array(x, list(i1, c(6:5, 6L), 3L, 2:1))
    checkIdentical(unname(a[i1, c(6:5, 6L), 3L, 2:1, drop=FALSE]),
                   unname(current))
    current <- extract_array(x, list(i1, c(6:5, 6L), integer(0), 2:1))
    checkIdentical(unname(a[i1, c(6:5, 6L), integer(0), 2:1]),
                   unname(current))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### test_* functions
###

test_DelayedSubset_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedSubset()
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(NULL), x@index)
    checkTrue(is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a)
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(NULL, 5:4))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, 5:4), x@index)
    checkIdentical(FALSE, is_noop(x))

    ## Test normalization of user-supplied 'Nindex' argument

    checkException(new_DelayedSubset(.TEST_MATRIX2a, 5:4),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX2a, list(5:4)),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX2a, list(NULL, 7:4)),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX2a, list(NULL, "zz")),
                   silent=silent)

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(0, c(5:4, 5.99, 6.99)))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(integer(0), c(5:4, 5:6)), x@index)
    checkIdentical(FALSE, is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(NULL, -3))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, c(1:2, 4:6)), x@index)
    checkIdentical(FALSE, is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(NULL, 1:6))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(-22, TRUE))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, NULL), x@index)
    checkTrue(is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(NULL, c("E", "B", "C", "C")))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(NULL, c(5L, 2L, 3L, 3L)), x@index)
    checkIdentical(FALSE, is_noop(x))

    x <- new_DelayedSubset(.TEST_MATRIX2a, list(IRanges(3:2, 5), IRanges(1, 6)))
    checkTrue(is(x, "DelayedSubset"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(c(3:5, 2:5), NULL), x@index)
    checkIdentical(FALSE, is_noop(x))

    checkException(new_DelayedSubset(.TEST_MATRIX2a, list(NULL, IRanges(0, 3))),
                   silent=silent)
    checkException(new_DelayedSubset(.TEST_MATRIX2a, list(NULL, IRanges(2, 7))),
                   silent=silent)
}

test_DelayedSubset_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedSubset(.TEST_MATRIX2a)

    .basic_checks_on_DelayedOp_with_DIM2(.TEST_MATRIX2a, x1)

    checkIdentical(FALSE, is_sparse(x1))

    ## 2. Ordinary array seed

    Nindex2 <- list(5:2, c(1:2, 5, 2:1))
    x2 <- new_DelayedSubset(.TEST_MATRIX2a, Nindex2)

    a2 <- .TEST_MATRIX2a[5:2, c(1:2, 5L, 2:1)]
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))

    current <- extract_array(x2, list(NULL, 4:2))
    checkIdentical(unname(a2[ , 4:2]), unname(current))
    current <- extract_array(x2, list(NULL, 4L))
    checkIdentical(unname(a2[ , 4L, drop=FALSE]), unname(current))
    current <- extract_array(x2, list(integer(0), NULL))
    checkIdentical(unname(a2[integer(0), ]), unname(current))

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedSubset(.TEST_SAS3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)
    .check_extract_sparse_array_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    ## 4. Sparse seed -- structural sparsity propagated

    i3 <- c(3:8, 1L)
    index4 <- list(i3, 6:5, NULL)
    x4 <- new_DelayedSubset(.TEST_SAS3, index4)

    a4 <- .TEST_ARRAY3[i3, 6:5, , drop=FALSE]
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(dimnames(a4), dimnames(x4))
    checkIdentical(a4, as.array(x4))

    current <- extract_array(x4, list(7:5, NULL, NULL))
    checkIdentical(unname(a4[7:5, , , drop=FALSE]), unname(current))
    current <- extract_array(x4, list(7:5, 2L, integer(0)))
    checkIdentical(unname(a4[7:5, 2L, integer(0), drop=FALSE]), unname(current))

    checkTrue(is_sparse(x4))
    sas4 <- extract_sparse_array(.TEST_SAS3, index4)
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

    i3 <- c(3:8, 1L)
    Nindex5 <- list(i3, c(6:5, 6L), NULL)  # Nindex5[[2]] contains duplicates!
    x5 <- new_DelayedSubset(.TEST_SAS3, Nindex5)

    a5 <- .TEST_ARRAY3[i3, c(6:5, 6L), , drop=FALSE]
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))

    checkIdentical(FALSE, is_sparse(x5))  # structural sparsity not propagated!
}

test_DelayedAperm_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedAperm()
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(1L, x@perm)
    checkTrue(is_noop(x))

    x <- new_DelayedAperm(.TEST_SAS3)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(1:3, x@perm)
    checkTrue(is_noop(x))

    x <- new_DelayedAperm(.TEST_SAS3, 3:1)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(3:1, x@perm)
    checkIdentical(FALSE, is_noop(x))

    checkException(new_DelayedAperm(.TEST_SAS3, "2"), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS3, integer(0)), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS3, NA_integer_), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS3, 9L), silent=silent)
    checkException(new_DelayedAperm(.TEST_SAS3, 1L), silent=silent)

    x <- new_DelayedAperm(.TEST_SAS3, 2:1)
    checkTrue(is(x, "DelayedAperm"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(2:1, x@perm)
    checkIdentical(FALSE, is_noop(x))
}

test_DelayedAperm_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedAperm(.TEST_MATRIX2a)

    .basic_checks_on_DelayedOp_with_DIM2(.TEST_MATRIX2a, x1)

    checkIdentical(FALSE, is_sparse(x1))

    ## 2. Ordinary array seed -- transposition

    x2 <- new_DelayedAperm(.TEST_MATRIX2a, 2:1)

    a2 <- t(.TEST_MATRIX2a)
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))

    current <- extract_array(x2, list(4:2, NULL))
    checkIdentical(unname(a2[4:2, ]), unname(current))
    current <- extract_array(x2, list(4L, NULL))
    checkIdentical(unname(a2[4L, , drop=FALSE]), unname(current))
    current <- extract_array(x2, list(NULL, integer(0)))
    checkIdentical(unname(a2[ , integer(0)]), unname(current))

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedAperm(.TEST_SAS3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)
    .check_extract_sparse_array_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    ## 4. Sparse seed -- transpose 1st and 3rd dims

    x4 <- new_DelayedAperm(.TEST_SAS3, 3:1)

    a4 <- aperm(.TEST_ARRAY3, 3:1)
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(dimnames(a4), dimnames(x4))
    checkIdentical(a4, as.array(x4))

    i3 <- c(3:8, 1L)
    current <- extract_array(x4, list(NULL, 6:5, i3))
    checkIdentical(unname(a4[ , 6:5, i3, drop=FALSE]), unname(current))
    current <- extract_array(x4, list(NULL, c(6:5, 6L), i3))
    checkIdentical(unname(a4[ , c(6:5, 6L), i3, drop=FALSE]), unname(current))
    current <- extract_array(x4, list(integer(0), c(6:5, 6L), i3))
    checkIdentical(unname(a4[integer(0), c(6:5, 6L), i3]), unname(current))

    checkTrue(is_sparse(x4))
    sas4 <- aperm(.TEST_SAS3, 3:1)
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x4, list(NULL, NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas4, current)
    current <- extract_sparse_array(x4, list(NULL, 6:5, i3))
    target <- extract_sparse_array(sas4, list(NULL, 6:5, i3))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x4, list(integer(0), 6:5, i3))
    target <- extract_sparse_array(sas4, list(integer(0), 6:5, i3))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)

    ## 5. Sparse seed -- transpose 1st and 2nd dims and drop 3rd dim

    x5 <- new_DelayedAperm(.TEST_SAS3, 2:1)

    a5 <- t(drop(.TEST_ARRAY3))
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))

    i3 <- c(3:8, 1L)
    current <- extract_array(x5, list(6:5, i3))
    checkIdentical(unname(a5[6:5, i3]), unname(current))
    current <- extract_array(x5, list(c(6:5, 6L), i3))
    checkIdentical(unname(a5[c(6:5, 6L), i3]), unname(current))
    current <- extract_array(x5, list(c(6:5, 6L), integer(0)))
    checkIdentical(unname(a5[c(6:5, 6L), integer(0)]), unname(current))

    checkTrue(is_sparse(x5))
    sas5 <- aperm(.TEST_SAS3, 2:1)
    ## The behavior of extract_sparse_array() is **undefined** when the
    ## subscripts in 'index' contain duplicates (see "extract_sparse_array()
    ## Terms of Use" in SparseArraySeed-class.R). So do NOT use such
    ## subscripts in the tests below.
    current <- extract_sparse_array(x5, list(NULL, NULL))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(sas5, current)
    current <- extract_sparse_array(x5, list(6:5, i3))
    target <- extract_sparse_array(sas5, list(6:5, i3))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
    current <- extract_sparse_array(x5, list(6:5, integer(0)))
    target <- extract_sparse_array(sas5, list(6:5, integer(0)))
    checkTrue(is(current, "SparseArraySeed"))
    checkIdentical(target, current)
}

test_DelayedUnaryIsoOpStack_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedUnaryIsoOpStack()
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(), x@OPS)
    #checkTrue(is_noop(x))  # no is_noop() yet for DelayedUnaryIsoOpStack objects

    x <- new_DelayedUnaryIsoOpStack(.TEST_SAS3)
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(list(), x@OPS)
    #checkTrue(is_noop(x))

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x <- new_DelayedUnaryIsoOpStack(.TEST_SAS3, OPS)
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(OPS, x@OPS)
    #checkIdentical(FALSE, is_noop(x))

    checkException(new_DelayedUnaryIsoOpStack(.TEST_SAS3, NULL),
                   silent=silent)
    checkException(new_DelayedUnaryIsoOpStack(.TEST_SAS3, list(NULL)),
                   silent=silent)
    checkException(new_DelayedUnaryIsoOpStack(.TEST_SAS3,
                                              list("not-an-existing-function")),
                   silent=silent)
}

test_DelayedUnaryIsoOpStack_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedUnaryIsoOpStack(.TEST_MATRIX2a)

    .basic_checks_on_DelayedOp_with_DIM2(.TEST_MATRIX2a, x1)

    checkIdentical(FALSE, is_sparse(x1))

    ## 2. Ordinary array seed -- 1 / (log(a)^2 + 1)

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x2 <- new_DelayedUnaryIsoOpStack(.TEST_MATRIX2a, OPS)

    a2 <- suppressWarnings(1 / (log(.TEST_MATRIX2a)^2 + 1))
    .basic_checks_on_DelayedOp_with_DIM2(a2, x2)

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedUnaryIsoOpStack(.TEST_SAS3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)
    .check_extract_sparse_array_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    ## 4. Sparse seed -- 1 / (log(a)^2 + 1)

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x4 <- new_DelayedUnaryIsoOpStack(.TEST_SAS3, OPS)

    a4 <- 1 / (log(.TEST_ARRAY3)^2 + 1)
    .basic_checks_on_DelayedOp_with_DIM3(a4, x4)
    .check_extract_sparse_array_on_DelayedOp_with_DIM3(a4, x4)

    ## 5. Sparse seed but structural sparsity NOT propagated because
    ##    the stack of operations doesn't preserve the zeros

    OPS <- list(function(a) cos(a),
                function(a) log(a^2 + 1))
    x5 <- new_DelayedUnaryIsoOpStack(.TEST_SAS3, OPS)

    a5 <- log(cos(.TEST_ARRAY3)^2 + 1)  # does not preserve the zeros
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))

    checkIdentical(FALSE, is_sparse(x5))  # structural sparsity not propagated!
}

test_DelayedUnaryIsoOpWithArgs_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedUnaryIsoOpWithArgs()
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(identity, x@OP)
    #checkTrue(is_noop(x))  # no is_noop() yet for DelayedUnaryIsoOpWithArgs
                           # objects

    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS3)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(identity, x@OP)
    #checkTrue(is_noop(x))

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    Rargs <- list(e2=e2)
    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS3, OP, Rargs=Rargs, Ralong=1L)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SAS3, seed(x))
    checkIdentical(OP, x@OP)
    checkIdentical(list(), x@Largs)
    checkIdentical(Rargs, x@Rargs)
    checkIdentical(integer(0), x@Lalong)
    checkIdentical(1L, x@Ralong)
    #checkIdentical(FALSE, is_noop(x))
}

test_DelayedUnaryIsoOpWithArgs_API <- function()
{
    ## Note that DelayedUnaryIsoOpWithArgs objects are NOT considered to
    ## propagate structural sparsity.

    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedUnaryIsoOpWithArgs(.TEST_MATRIX2a)

    .basic_checks_on_DelayedOp_with_DIM2(.TEST_MATRIX2a, x1)

    checkIdentical(FALSE, is_sparse(x1))

    ## 2. Ordinary array seed

    OP <- `/`
    e2 <- rowSums(.TEST_MATRIX2a)
    x2 <- new_DelayedUnaryIsoOpWithArgs(.TEST_MATRIX2a, OP,
                                        Rargs=list(e2=e2), Ralong=1L)

    a2 <- .TEST_MATRIX2a / e2
    .basic_checks_on_DelayedOp_with_DIM2(a2, x2)

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    checkIdentical(TRUE, is_sparse(x3))

    ## 4. Sparse seed

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    x4 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS3, OP,
                                        Rargs=list(e2=e2), Ralong=1L)

    a4 <- .TEST_ARRAY3 <= e2
    .basic_checks_on_DelayedOp_with_DIM3(a4, x4)

    checkIdentical(FALSE, is_sparse(x4))
}

test_DelayedSetDimnames_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedSetDimnames()
    checkTrue(is(x, "DelayedSetDimnames"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(.INHERIT_FROM_SEED), x@dimnames)
    checkTrue(is_noop(x))

    x <- new_DelayedSetDimnames(.TEST_MATRIX2a)
    checkTrue(is(x, "DelayedSetDimnames"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(rep(list(.INHERIT_FROM_SEED), 2), x@dimnames)
    checkTrue(is_noop(x))

    x <- new_DelayedSetDimnames(.TEST_MATRIX2a, dimnames(.TEST_MATRIX2a))
    checkTrue(is(x, "DelayedSetDimnames"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(rep(list(.INHERIT_FROM_SEED), 2), x@dimnames)
    checkTrue(is_noop(x))

    x <- new_DelayedSetDimnames(.TEST_MATRIX2a, NULL)
    checkTrue(is(x, "DelayedSetDimnames"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(.INHERIT_FROM_SEED, NULL), x@dimnames)
    checkIdentical(FALSE, is_noop(x))
    x2 <- new_DelayedSetDimnames(.TEST_MATRIX2a, list(NULL))
    checkIdentical(x, x2)
    x3 <- new_DelayedSetDimnames(.TEST_MATRIX2a, list(NULL, NULL))
    checkIdentical(x, x3)

    x <- new_DelayedSetDimnames(.TEST_MATRIX2a, list(letters[1:5], NULL))
    checkTrue(is(x, "DelayedSetDimnames"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_MATRIX2a, seed(x))
    checkIdentical(list(letters[1:5], NULL), x@dimnames)
    checkIdentical(FALSE, is_noop(x))
    x2 <- new_DelayedSetDimnames(.TEST_MATRIX2a, list(letters[1:5]))
    checkIdentical(x, x2)

    checkException(new_DelayedSetDimnames(.TEST_MATRIX2a, letters),
                   silent=silent)
    checkException(new_DelayedSetDimnames(.TEST_MATRIX2a, list(NULL, NULL, NULL)),
                   silent=silent)
    checkException(new_DelayedSetDimnames(.TEST_MATRIX2a, list(ls, NULL)),
                   silent=silent)
    checkException(new_DelayedSetDimnames(.TEST_MATRIX2a, list(letters, NULL)),
                   silent=silent)
}

test_DelayedSetDimnames_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedSetDimnames(.TEST_MATRIX2a)

    .basic_checks_on_DelayedOp_with_DIM2(.TEST_MATRIX2a, x1)

    checkIdentical(FALSE, is_sparse(x1))

    ## 2. Ordinary array seed

    dimnames2 <- list(letters[1:5], NULL)
    x2 <- new_DelayedSetDimnames(.TEST_MATRIX2a, dimnames2)

    a2 <- .TEST_MATRIX2a
    dimnames(a2) <- dimnames2
    .basic_checks_on_DelayedOp_with_DIM2(a2, x2)

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedSetDimnames(.TEST_SAS3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    checkTrue(is_sparse(x3))

    ## 4. Sparse seed

    dimnames4 <- list(NULL, sprintf("ID%03d", seq_len(ncol(.TEST_SAS3))), NULL)
    x4 <- new_DelayedSetDimnames(.TEST_SAS3, dimnames4)

    a4 <- .TEST_ARRAY3
    dimnames(a4) <- dimnames4
    .basic_checks_on_DelayedOp_with_DIM3(a4, x4)

    checkTrue(is_sparse(x4))
}

test_DelayedNaryIsoOp_constructor <- function(silent=FALSE)
{
    ## We also test nseed(), seed(), and is_noop()

    x <- new_DelayedNaryIsoOp()
    checkTrue(is(x, "DelayedNaryIsoOp"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkException(seed(x), silent=silent)
    checkIdentical(identity, x@OP)
    checkIdentical(list(new("array")), x@seeds)
    checkIdentical(list(), x@Rargs)
    checkException(is_noop(x), silent=silent)

    x <- new_DelayedNaryIsoOp("/", .TEST_ARRAY4a, .TEST_SAS4c)
    checkTrue(is(x, "DelayedNaryIsoOp"))
    checkTrue(validObject(x))
    checkIdentical(2L, nseed(x))
    checkException(seed(x), silent=silent)
    checkIdentical(`/`, x@OP)
    checkIdentical(list(.TEST_ARRAY4a, .TEST_SAS4c), x@seeds)
    checkIdentical(list(), x@Rargs)
    checkException(is_noop(x), silent=silent)

    OP <- function(a1, a2, a3) (a1 + a2 + a3) / 3  # ternary mean
    x <- new_DelayedNaryIsoOp(OP, .TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c)
    checkTrue(is(x, "DelayedNaryIsoOp"))
    checkTrue(validObject(x))
    checkIdentical(3L, nseed(x))
    checkException(seed(x), silent=silent)
    checkIdentical(OP, x@OP)
    checkIdentical(list(.TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c), x@seeds)
    checkIdentical(list(), x@Rargs)
    checkException(is_noop(x), silent=silent)

    ## Two alternate ways to represent the above DelayedNaryIsoOp.

    OP <- function(a1, a2, a3, d=3) (a1 + a2 + a3) / d 
    x <- new_DelayedNaryIsoOp(OP, .TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c)
    checkTrue(is(x, "DelayedNaryIsoOp"))
    checkTrue(validObject(x))
    checkIdentical(3L, nseed(x))
    checkException(seed(x), silent=silent)
    checkIdentical(OP, x@OP)
    checkIdentical(list(.TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c), x@seeds)
    checkIdentical(list(), x@Rargs)
    checkException(is_noop(x), silent=silent)

    OP <- function(a1, a2, a3, d) (a1 + a2 + a3) / d 
    Rargs <- list(d=3)
    x <- new_DelayedNaryIsoOp(OP, .TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c,
                                  Rargs=Rargs)
    checkTrue(is(x, "DelayedNaryIsoOp"))
    checkTrue(validObject(x))
    checkIdentical(3L, nseed(x))
    checkException(seed(x), silent=silent)
    checkIdentical(OP, x@OP)
    checkIdentical(list(.TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c), x@seeds)
    checkIdentical(Rargs, x@Rargs)
    checkException(is_noop(x), silent=silent)

    checkException(new_DelayedNaryIsoOp(NULL),
                   silent=silent)
    checkException(new_DelayedNaryIsoOp(list(NULL)),
                   silent=silent)
    checkException(new_DelayedNaryIsoOp("not-an-existing-function"),
                   silent=silent)
    checkException(new_DelayedNaryIsoOp("<=", array(dim=4:2), array(dim=2:4)),
                   silent=silent)
}

test_DelayedNaryIsoOp_API <- function()
{
    ## 1. Unary op

    OP1 <- function(a) 1 / (log(a)^2 + 1)
    x1 <- new_DelayedNaryIsoOp(OP1, .TEST_SAS4c)

    a1 <- OP1(.TEST_ARRAY4c)
    .basic_checks_on_DelayedOp_with_DIM4(a1, x1)

    checkTrue(is_sparse(x1))

    ## 2. Binary op

    x2 <- new_DelayedNaryIsoOp(`/`, .TEST_ARRAY4a, .TEST_SAS4c)

    a2 <- .TEST_ARRAY4a / .TEST_ARRAY4c
    .basic_checks_on_DelayedOp_with_DIM4(a2, x2)

    checkIdentical(FALSE, is_sparse(x2))

    ## 3. Ternary op with right args

    ## Ternary array sum with weights on the 2nd and 3rd arrays (weights
    ## should be single values)
    OP3 <- function(a1, a2, a3, w2, w3) a1 + w2*a2 + w3*a3
    Rargs3 <- list(w2=0.5, w3=100L)
    x3 <- new_DelayedNaryIsoOp(OP3, .TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_SAS4c,
                                    Rargs=Rargs3)

    a3 <- OP3(.TEST_ARRAY4a, .TEST_ARRAY4b, .TEST_ARRAY4c, w2=0.5, w3=100L)
    .basic_checks_on_DelayedOp_with_DIM4(a3, x3)

    checkIdentical(FALSE, is_sparse(x3))

    ## 4-11. "Ops" group generics with propagation of structural sparsity

    x4 <- new_DelayedNaryIsoOp("+", .TEST_SAS4c, .TEST_SAS4d)

    a4 <- .TEST_ARRAY4c + .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a4, x4)

    checkTrue(is_sparse(x4))

    x5 <- new_DelayedNaryIsoOp("-", .TEST_SAS4c, .TEST_SAS4d)

    a5 <- .TEST_ARRAY4c - .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a5, x5)

    checkTrue(is_sparse(x5))

    x6 <- new_DelayedNaryIsoOp("*", .TEST_SAS4c, .TEST_SAS4d)

    a6 <- .TEST_ARRAY4c * .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a6, x6)

    checkTrue(is_sparse(x6))

    x7 <- new_DelayedNaryIsoOp(">", .TEST_SAS4c, .TEST_SAS4d)

    a7 <- .TEST_ARRAY4c > .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a7, x7)

    checkTrue(is_sparse(x7))

    x8 <- new_DelayedNaryIsoOp("<", .TEST_SAS4c, .TEST_SAS4d)

    a8 <- .TEST_ARRAY4c < .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a8, x8)

    checkTrue(is_sparse(x8))

    x9 <- new_DelayedNaryIsoOp("!=", .TEST_SAS4c, .TEST_SAS4d)

    a9 <- .TEST_ARRAY4c != .TEST_ARRAY4d
    .basic_checks_on_DelayedOp_with_DIM4(a9, x9)

    checkTrue(is_sparse(x9))

    x10 <- new_DelayedNaryIsoOp("&", .TEST_SAS2b, .TEST_SAS2c)

    a10 <- .TEST_ARRAY2b & .TEST_ARRAY2c
    .basic_checks_on_DelayedOp_with_DIM2(a10, x10)

    checkTrue(is_sparse(x10))

    x11 <- new_DelayedNaryIsoOp("|", .TEST_SAS2b, .TEST_SAS2c)

    a11 <- .TEST_ARRAY2b | .TEST_ARRAY2c
    .basic_checks_on_DelayedOp_with_DIM2(a11, x11)

    checkTrue(is_sparse(x11))

    ## 12. "pmax2" with propagation of structural sparsity

    # TODO

    ## 13. "pmin2" with propagation of structural sparsity

    # TODO
}

