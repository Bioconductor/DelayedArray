new_DelayedSubset <- DelayedArray:::new_DelayedSubset
new_DelayedAperm <- DelayedArray:::new_DelayedAperm
new_DelayedUnaryIsoOpStack <- DelayedArray:::new_DelayedUnaryIsoOpStack
new_DelayedUnaryIsoOpWithArgs <- DelayedArray:::new_DelayedUnaryIsoOpWithArgs
new_DelayedDimnames <- DelayedArray:::new_DelayedDimnames
.INHERIT_FROM_SEED <- DelayedArray:::.INHERIT_FROM_SEED

.TEST_MATRIX1 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6,
                        dimnames=list(NULL, LETTERS[1:6]))

.TEST_SAS2 <- SparseArraySeed(dim=c(40, 100, 1),
                              aind=arrayInd(1:50, c(10, 5, 1)),
                              nzdata=1:50)
.TEST_ARRAY2 <- sparse2dense(.TEST_SAS2)

.basic_checks_on_DelayedOp_with_MATRIX1_dim <- function(a, x)
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

.basic_checks_on_DelayedOp_with_ARRAY2_dim <- function(a, x)
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

.check_extract_sparse_array_on_DelayedOp_with_ARRAY2_dim <- function(a, x)
{
    checkTrue(isSparse(x))

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### test_* functions
###

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
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedSubset(.TEST_MATRIX1)

    .basic_checks_on_DelayedOp_with_MATRIX1_dim(.TEST_MATRIX1, x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed

    Nindex2 <- list(5:2, c(1:2, 5, 2:1))
    x2 <- new_DelayedSubset(.TEST_MATRIX1, Nindex2)

    a2 <- .TEST_MATRIX1[5:2, c(1:2, 5L, 2:1)]
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))

    current <- extract_array(x2, list(NULL, 4:2))
    checkIdentical(unname(a2[ , 4:2]), unname(current))
    current <- extract_array(x2, list(NULL, 4L))
    checkIdentical(unname(a2[ , 4L, drop=FALSE]), unname(current))
    current <- extract_array(x2, list(integer(0), NULL))
    checkIdentical(unname(a2[integer(0), ]), unname(current))

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedSubset(.TEST_SAS2)

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)
    .check_extract_sparse_array_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)

    ## 4. Sparse seed -- structural sparsity propagated

    i3 <- c(3:8, 1L)
    index4 <- list(i3, 6:5, NULL)
    x4 <- new_DelayedSubset(.TEST_SAS2, index4)

    a4 <- .TEST_ARRAY2[i3, 6:5, , drop=FALSE]
    checkIdentical(dim(a4), dim(x4))
    checkIdentical(dimnames(a4), dimnames(x4))
    checkIdentical(a4, as.array(x4))

    current <- extract_array(x4, list(7:5, NULL, NULL))
    checkIdentical(unname(a4[7:5, , , drop=FALSE]), unname(current))
    current <- extract_array(x4, list(7:5, 2L, integer(0)))
    checkIdentical(unname(a4[7:5, 2L, integer(0), drop=FALSE]), unname(current))

    checkTrue(isSparse(x4))
    sas4 <- extract_sparse_array(.TEST_SAS2, index4)
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
    x5 <- new_DelayedSubset(.TEST_SAS2, Nindex5)

    a5 <- .TEST_ARRAY2[i3, c(6:5, 6L), , drop=FALSE]
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
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedAperm(.TEST_MATRIX1)

    .basic_checks_on_DelayedOp_with_MATRIX1_dim(.TEST_MATRIX1, x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed -- transposition

    x2 <- new_DelayedAperm(.TEST_MATRIX1, 2:1)

    a2 <- t(.TEST_MATRIX1)
    checkIdentical(dim(a2), dim(x2))
    checkIdentical(dimnames(a2), dimnames(x2))
    checkIdentical(a2, as.array(x2))

    current <- extract_array(x2, list(4:2, NULL))
    checkIdentical(unname(a2[4:2, ]), unname(current))
    current <- extract_array(x2, list(4L, NULL))
    checkIdentical(unname(a2[4L, , drop=FALSE]), unname(current))
    current <- extract_array(x2, list(NULL, integer(0)))
    checkIdentical(unname(a2[ , integer(0)]), unname(current))

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedAperm(.TEST_SAS2)

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)
    .check_extract_sparse_array_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)

    ## 4. Sparse seed -- transpose 1st and 3rd dims

    x4 <- new_DelayedAperm(.TEST_SAS2, 3:1)

    a4 <- aperm(.TEST_ARRAY2, 3:1)
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

    checkTrue(isSparse(x4))
    sas4 <- aperm(.TEST_SAS2, 3:1)
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

    x5 <- new_DelayedAperm(.TEST_SAS2, 2:1)

    a5 <- t(drop(.TEST_ARRAY2))
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

    checkTrue(isSparse(x5))
    sas5 <- aperm(.TEST_SAS2, 2:1)
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
    ## We also test seed() and isNoOp()

    x <- new_DelayedUnaryIsoOpStack()
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(), x@OPS)
    #checkTrue(isNoOp(x))  # no isNoOp() yet for DelayedUnaryIsoOpStack objects

    x <- new_DelayedUnaryIsoOpStack(.TEST_SAS2)
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(list(), x@OPS)
    #checkTrue(isNoOp(x))

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x <- new_DelayedUnaryIsoOpStack(.TEST_SAS2, OPS)
    checkTrue(is(x, "DelayedUnaryIsoOpStack"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(OPS, x@OPS)
    #checkIdentical(FALSE, isNoOp(x))
}

test_DelayedUnaryIsoOpStack_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedUnaryIsoOpStack(.TEST_MATRIX1)

    .basic_checks_on_DelayedOp_with_MATRIX1_dim(.TEST_MATRIX1, x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed -- 1 / (log(a)^2 + 1)

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x2 <- new_DelayedUnaryIsoOpStack(.TEST_MATRIX1, OPS)

    a2 <- suppressWarnings(1 / (log(.TEST_MATRIX1)^2 + 1))
    .basic_checks_on_DelayedOp_with_MATRIX1_dim(a2, x2)

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedUnaryIsoOpStack(.TEST_SAS2)

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)
    .check_extract_sparse_array_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)

    ## 4. Sparse seed -- 1 / (log(a)^2 + 1)

    OPS <- list(function(a) log(a),
                function(a) a^2 + 1,
                function(a) 1 / a)
    x4 <- new_DelayedUnaryIsoOpStack(.TEST_SAS2, OPS)

    a4 <- 1 / (log(.TEST_ARRAY2)^2 + 1)
    .basic_checks_on_DelayedOp_with_ARRAY2_dim(a4, x4)
    .check_extract_sparse_array_on_DelayedOp_with_ARRAY2_dim(a4, x4)

    ## 5. Sparse seed but structural sparsity NOT propagated because
    ##    the stack of operations doesn't preserve the zeroes

    OPS <- list(function(a) cos(a),
                function(a) log(a^2 + 1))
    x5 <- new_DelayedUnaryIsoOpStack(.TEST_SAS2, OPS)

    a5 <- log(cos(.TEST_ARRAY2)^2 + 1)  # does not preserve the zeroes
    checkIdentical(dim(a5), dim(x5))
    checkIdentical(dimnames(a5), dimnames(x5))
    checkIdentical(a5, as.array(x5))

    checkIdentical(FALSE, isSparse(x5))  # structural sparsity not propagated!
}

test_DelayedUnaryIsoOpWithArgs_constructor <- function(silent=FALSE)
{
    ## We also test seed() and isNoOp()

    x <- new_DelayedUnaryIsoOpWithArgs()
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(identity, x@OP)
    #checkTrue(isNoOp(x))  # no isNoOp() yet for DelayedUnaryIsoOpWithArgs
                           # objects

    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS2)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(identity, x@OP)
    #checkTrue(isNoOp(x))

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    Rargs <- list(e2=e2)
    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS2, OP, Rargs=Rargs, Ralong=1L)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_SAS2, seed(x))
    checkIdentical(OP, x@OP)
    checkIdentical(list(), x@Largs)
    checkIdentical(Rargs, x@Rargs)
    checkIdentical(integer(0), x@Lalong)
    checkIdentical(1L, x@Ralong)
    #checkIdentical(FALSE, isNoOp(x))
}

test_DelayedUnaryIsoOpWithArgs_API <- function()
{
    ## Note that DelayedUnaryIsoOpWithArgs objects are NOT considered to
    ## propagate structural sparsity.

    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedUnaryIsoOpWithArgs(.TEST_MATRIX1)

    .basic_checks_on_DelayedOp_with_MATRIX1_dim(.TEST_MATRIX1, x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed

    OP <- `/`
    e2 <- rowSums(.TEST_MATRIX1)
    x2 <- new_DelayedUnaryIsoOpWithArgs(.TEST_MATRIX1, OP,
                                        Rargs=list(e2=e2), Ralong=1L)

    a2 <- .TEST_MATRIX1 / e2
    .basic_checks_on_DelayedOp_with_MATRIX1_dim(a2, x2)

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS2)

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)

    checkIdentical(FALSE, isSparse(x3))

    ## 4. Sparse seed

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    x4 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SAS2, OP,
                                        Rargs=list(e2=e2), Ralong=1L)

    a4 <- .TEST_ARRAY2 <= e2

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(a4, x4)

    checkIdentical(FALSE, isSparse(x4))
}

test_DelayedDimnames_constructor <- function(silent=FALSE)
{
    ## We also test seed() and isNoOp()

    x <- new_DelayedDimnames()
    checkTrue(is(x, "DelayedDimnames"))
    checkTrue(validObject(x))
    checkIdentical(new("array"), seed(x))
    checkIdentical(list(.INHERIT_FROM_SEED), x@dimnames)
    checkTrue(isNoOp(x))

    x <- new_DelayedDimnames(.TEST_MATRIX1)
    checkTrue(is(x, "DelayedDimnames"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(rep(list(.INHERIT_FROM_SEED), 2), x@dimnames)
    checkTrue(isNoOp(x))

    x <- new_DelayedDimnames(.TEST_MATRIX1, dimnames(.TEST_MATRIX1))
    checkTrue(is(x, "DelayedDimnames"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(rep(list(.INHERIT_FROM_SEED), 2), x@dimnames)
    checkTrue(isNoOp(x))

    x <- new_DelayedDimnames(.TEST_MATRIX1, NULL)
    checkTrue(is(x, "DelayedDimnames"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(.INHERIT_FROM_SEED, NULL), x@dimnames)
    checkIdentical(FALSE, isNoOp(x))

    x <- new_DelayedDimnames(.TEST_MATRIX1, list(letters[1:5], NULL))
    checkTrue(is(x, "DelayedDimnames"))
    checkTrue(validObject(x))
    checkIdentical(.TEST_MATRIX1, seed(x))
    checkIdentical(list(letters[1:5], NULL), x@dimnames)
    checkIdentical(FALSE, isNoOp(x))

    checkException(new_DelayedDimnames(.TEST_MATRIX1, letters),
                   silent=silent)
    checkException(new_DelayedDimnames(.TEST_MATRIX1, list(NULL)),
                   silent=silent)
    checkException(new_DelayedDimnames(.TEST_MATRIX1, list(NULL, NULL, NULL)),
                   silent=silent)
    checkException(new_DelayedDimnames(.TEST_MATRIX1, list(1:5, NULL)),
                   silent=silent)
    checkException(new_DelayedDimnames(.TEST_MATRIX1, list(letters, NULL)),
                   silent=silent)
}

test_DelayedDimnames_API <- function()
{
    ## 1. Ordinary array seed -- no-op

    x1 <- new_DelayedDimnames(.TEST_MATRIX1)

    .basic_checks_on_DelayedOp_with_MATRIX1_dim(.TEST_MATRIX1, x1)

    checkIdentical(FALSE, isSparse(x1))

    ## 2. Ordinary array seed

    dimnames2 <- list(letters[1:5], NULL)
    x2 <- new_DelayedDimnames(.TEST_MATRIX1, dimnames2)

    a2 <- .TEST_MATRIX1
    dimnames(a2) <- dimnames2
    .basic_checks_on_DelayedOp_with_MATRIX1_dim(a2, x2)

    checkIdentical(FALSE, isSparse(x2))

    ## 3. Sparse seed -- no-op

    x3 <- new_DelayedDimnames(.TEST_SAS2)

    .basic_checks_on_DelayedOp_with_ARRAY2_dim(.TEST_ARRAY2, x3)

    checkTrue(isSparse(x3))

    ## 4. Sparse seed

    dimnames4 <- list(NULL, sprintf("ID%03d", seq_len(ncol(.TEST_SAS2))), NULL)
    x4 <- new_DelayedDimnames(.TEST_SAS2, dimnames4)

    a4 <- .TEST_ARRAY2
    dimnames(a4) <- dimnames4
    .basic_checks_on_DelayedOp_with_ARRAY2_dim(a4, x4)

    checkTrue(isSparse(x4))
}

