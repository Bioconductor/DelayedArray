#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

.TEST_matrices <- list(
    matrix(1:15, nrow=3, ncol=5,
           dimnames=list(NULL, paste0("M1y", 1:5))),
    matrix(101:135, nrow=7, ncol=5,
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5))),
    matrix(1001:1025, nrow=5, ncol=5,
           dimnames=list(paste0("M3x", 1:5), NULL))
)

.TEST_arrays <- list(
    array(1:60, c(3, 5, 4),
           dimnames=list(NULL, paste0("M1y", 1:5), NULL)),
    array(101:240, c(7, 5, 4),
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL)),
    array(10001:10100, c(5, 5, 4),
           dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))
)

test_DelayedMatrix_rbind_cbind <- function()
{
    m1 <- .TEST_matrices[[1]]
    m2 <- .TEST_matrices[[2]]
    m3 <- .TEST_matrices[[3]]
    M1 <- realize(m1)
    M2 <- realize(m2)
    M3 <- realize(m3)

    target <- rbind(a=m1, b=m2, c=m3)
    current <- rbind(a=M1, b=M2, c=M3)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    current <- cbind(a=t(M1), b=t(M2), c=t(M3))
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(t(target), as.matrix(current))

    ## unary form

    target <- rbind(a=m1)
    current <- rbind(a=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- cbind(a=m1)
    current <- cbind(a=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    ## with empty matrices

    m1 <- matrix(nrow=0, ncol=3, dimnames=list(NULL, letters[1:3]))
    m2 <- matrix(1:15, ncol=3, dimnames=list(NULL, LETTERS[1:3]))
    M1 <- realize(m1)
    M2 <- realize(m2)

    target <- rbind(a=m1, a=m2)
    current <- rbind(a=M1, b=M2)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- rbind(a=m2, a=m1)
    current <- rbind(a=M2, b=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- rbind(a=m1, a=m1)
    current <- rbind(a=M1, b=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))
}

test_DelayedArray_arbind <- function()
{
    TEST_hdf5arrays <- lapply(.TEST_arrays, realize)

    target <- do.call(arbind, .TEST_arrays)
    current <- do.call(arbind, TEST_hdf5arrays)
    checkTrue(is(current, "DelayedArray"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.array(current))

    ## For some mysterious reason, the code below fails in the context
    ## of 'R CMD check' but not when running the tests interactively with
    ## DelayedArray:::.test().
    #check_2D_slice <- function(k) {
    #    slices <- lapply(lapply(TEST_hdf5arrays, `[`, , , k), drop)
    #    target_slice <- do.call(rbind, slices)
    #    checkIdentical(as.matrix(target_slice), as.matrix(current[ , , k]))
    #}
    #for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)
}

test_DelayedArray_acbind <- function()
{
    ## transpose the 1st 2 dimensions
    arrays <- lapply(.TEST_arrays,
        function(a) {
            a_dimnames <- dimnames(a)
            dim(a)[1:2] <- dim(a)[2:1]
            a_dimnames[1:2] <- a_dimnames[2:1]
            dimnames(a) <- a_dimnames
            a
    })
    TEST_hdf5arrays <- lapply(arrays, realize)

    target <- do.call(acbind, arrays)
    current <- do.call(acbind, TEST_hdf5arrays)
    checkTrue(is(current, "DelayedArray"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.array(current))

    ## For some mysterious reason, the code below fails in the context
    ## of 'R CMD check' but not when running the tests interactively with
    ## DelayedArray:::.test().
    #check_2D_slice <- function(k) {
    #    slices <- lapply(lapply(TEST_hdf5arrays, `[`, , , k), drop)
    #    target_slice <- do.call(cbind, slices)
    #    checkIdentical(as.matrix(target_slice), as.matrix(current[ , , k]))
    #}
    #for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)
}

