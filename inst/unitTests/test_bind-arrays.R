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

test_arbind_on_matrices <- function()
{
    ## on non-empty matrices
    matrices <- .TEST_matrices
    target <- do.call(rbind, matrices)
    current <- do.call(arbind, matrices)
    checkIdentical(target, current)

    ## on empty matrices
    m1 <- matrix(nrow=0, ncol=3, dimnames=list(NULL, letters[1:3]))
    m2 <- matrix(1:15, ncol=3, dimnames=list(NULL, LETTERS[1:3]))

    target <- do.call(rbind, list(m1, m2))
    current <- do.call(arbind, list(m1, m2))
    checkIdentical(target, current)

    target <- do.call(rbind, list(m2, m1))
    current <- do.call(arbind, list(m2, m1))
    checkIdentical(target, current)

    target <- do.call(rbind, list(m1, m1))
    current <- do.call(arbind, list(m1, m1))
    checkIdentical(target, current)

    ## on matrices of type "list"
    m3 <- matrix(list(), nrow=0, ncol=5)
    m4 <- matrix(list(10:9, NULL, letters[1:3], TRUE, raw()), nrow=4, ncol=5)

    target <- do.call(rbind, list(m3, m4))
    current <- do.call(arbind, list(m3, m4))
    checkIdentical(target, current)

    target <- do.call(rbind, list(m4, m3))
    current <- do.call(arbind, list(m4, m3))
    checkIdentical(target, current)

    target <- do.call(rbind, list(m3, m3))
    current <- do.call(arbind, list(m3, m3))
    checkIdentical(target, current)

    matrices <- c(matrices, list(m3), list(m4))
    target <- do.call(rbind, matrices)
    current <- do.call(arbind, matrices)
    checkIdentical(target, current)
}

test_arbind_on_arrays <- function()
{
    ## on arbitrary arrays
    current <- do.call(arbind, .TEST_arrays)
    check_2D_slice <- function(k) {
        slices <- lapply(.TEST_arrays, `[`, , , k)
        target_slice <- do.call(rbind, slices)
        checkIdentical(target_slice, current[ , , k])
    }
    for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)

    ## on 1D arrays
    a1 <- array(11:15, 5, dimnames=list(LETTERS[1:5]))
    checkIdentical(a1, arbind(a1))                      # unary
    b1 <- array(letters[1:3])
    target <- array(c(a1, b1), 8, dimnames=list(c(LETTERS[1:5], rep("", 3))))
    checkIdentical(target, arbind(a1, b1))              # binary
    target <- array(c(b1, a1), 8, dimnames=list(c(rep("", 3), LETTERS[1:5])))
    checkIdentical(target, arbind(b1, a1))              # binary
    a1b1a1 <- arbind(a1, b1, a1)                        # ternary
    checkIdentical(arbind(a1, arbind(b1, a1)), a1b1a1)
    checkIdentical(arbind(arbind(a1, b1), a1), a1b1a1)
}

test_acbind_on_matrices <- function()
{
    ## on non-empty matrices
    matrices <- lapply(.TEST_matrices, t)
    target <- do.call(cbind, matrices)
    current <- do.call(acbind, matrices)
    checkIdentical(target, current)

    ## on empty matrices
    m1 <- matrix(nrow=3, ncol=0, dimnames=list(letters[1:3], NULL))
    m2 <- matrix(1:15, nrow=3, dimnames=list(LETTERS[1:3], NULL))

    target <- do.call(cbind, list(m1, m2))
    current <- do.call(acbind, list(m1, m2))
    checkIdentical(target, current)

    target <- do.call(cbind, list(m2, m1))
    current <- do.call(acbind, list(m2, m1))
    checkIdentical(target, current)

    target <- do.call(cbind, list(m1, m1))
    current <- do.call(acbind, list(m1, m1))
    checkIdentical(target, current)

    ## on matrices of type "list"
    m3 <- matrix(list(), nrow=5, ncol=0)
    m4 <- matrix(list(), nrow=5, ncol=4)

    target <- do.call(cbind, list(m3, m4))
    current <- do.call(acbind, list(m3, m4))
    checkIdentical(target, current)

    target <- do.call(cbind, list(m4, m3))
    current <- do.call(acbind, list(m4, m3))
    checkIdentical(target, current)

    target <- do.call(cbind, list(m3, m3))
    current <- do.call(acbind, list(m3, m3))
    checkIdentical(target, current)

    matrices <- c(matrices, list(m3), list(m4))
    target <- do.call(cbind, matrices)
    current <- do.call(acbind, matrices)
    checkIdentical(target, current)
}

test_acbind_on_arrays <- function()
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

    current <- do.call(acbind, arrays)
    check_2D_slice <- function(k) {
        slices <- lapply(arrays, `[`, , , k)
        target_slice <- do.call(cbind, slices)
        checkIdentical(target_slice, current[ , , k])
    }
    for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)

    ## acbind() is not supported on 1D arrays
    a1 <- array(11:15, 5, dimnames=list(LETTERS[1:5]))
    checkException(acbind(a1))          # unary
    b1 <- array(letters[1:3])
    checkException(acbind(a1, b1))      # binary
    checkException(acbind(b1, a1))      # binary
    checkException(acbind(a1, b1, a1))  # ternary
}

