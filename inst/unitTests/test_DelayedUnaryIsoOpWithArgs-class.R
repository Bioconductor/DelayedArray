new_DelayedUnaryIsoOpWithArgs <- DelayedArray:::new_DelayedUnaryIsoOpWithArgs

.TEST_MATRIX2a <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6,
                         dimnames=list(NULL, LETTERS[1:6]))

.TEST_SVT3 <- SVT_SparseArray(dim=c(40, 100, 1))
.TEST_SVT3[Lindex2Mindex(1:50, c(10, 5, 1))] <- 1:50

.TEST_ARRAY3 <- as.array(.TEST_SVT3)


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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### test_* functions
###

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

    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SVT3)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SVT3, seed(x))
    checkIdentical(identity, x@OP)
    #checkTrue(is_noop(x))

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    Rargs <- list(e2=e2)
    x <- new_DelayedUnaryIsoOpWithArgs(.TEST_SVT3, OP, Rargs=Rargs, Ralong=1L)
    checkTrue(is(x, "DelayedUnaryIsoOpWithArgs"))
    checkTrue(validObject(x))
    checkIdentical(1L, nseed(x))
    checkIdentical(.TEST_SVT3, seed(x))
    checkIdentical(OP, x@OP)
    checkIdentical(list(), x@Largs)
    checkIdentical(Rargs, x@Rargs)
    checkIdentical(integer(0), x@Lalong)
    checkIdentical(1L, x@Ralong)
    #checkIdentical(FALSE, is_noop(x))
}

test_DelayedUnaryIsoOpWithArgs_API <- function()
{
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

    x3 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SVT3)

    .basic_checks_on_DelayedOp_with_DIM3(.TEST_ARRAY3, x3)

    checkIdentical(TRUE, is_sparse(x3))

    ## 4. Sparse seed

    OP <- `<=`
    e2 <- rep(c(5, 10), 20)
    x4 <- new_DelayedUnaryIsoOpWithArgs(.TEST_SVT3, OP,
                                        Rargs=list(e2=e2), Ralong=1L)

    a4 <- .TEST_ARRAY3 <= e2
    .basic_checks_on_DelayedOp_with_DIM3(a4, x4)

    checkIdentical(FALSE, is_sparse(x4))

    ## 5. Sparse seed with a vector-like argument

    x5 <- new_DelayedUnaryIsoOpWithArgs(
        ConstantArraySeed(c(2, 4), value=0), `+`,
        Rargs=list(e2=c(0, 0)), Ralong=1L
    )
    checkTrue(is_sparse(x5))
    a5 <- matrix(0, nrow=2, ncol=4)
    checkIdentical(a5, as.array(x5))
    checkIdentical(a5, as.array(as(x5, "SVT_SparseArray")))
}
