
test_native_cbind <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)
    Y <- Matrix::rsparsematrix(1000, 50, 0.1)
    Dy <- DelayedArray(Y)

    Dr <- cbind(Dx, Dy)
    ref <- as.matrix(cbind(X, Y))
    dimnames(ref) <- NULL
    checkIdentical(as.matrix(Dr), ref) 

    # Works with a mix of matrix types.
    Dy2 <- DelayedArray(as.matrix(Y))
    Dr2 <- cbind(Dx, Dy2)
    checkIdentical(as.matrix(Dr2), ref) 

    # Works with >2 matrices.
    Dr3 <- cbind(Dx, Dy, Dx, Dy)
    ref <- as.matrix(cbind(X, Y, X, Y))
    dimnames(ref) <- NULL
    checkIdentical(as.matrix(Dr3), ref) 
}

test_native_rbind <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)
    Y <- Matrix::rsparsematrix(5000, 100, 0.1)
    Dy <- DelayedArray(Y)

    Dr <- rbind(Dx, Dy)
    ref <- as.matrix(rbind(X, Y))
    dimnames(ref) <- NULL
    checkIdentical(as.matrix(Dr), ref) 

    # Works with a mix of matrix types.
    Dy2 <- DelayedArray(as.matrix(Y))
    Dr2 <- rbind(Dx, Dy2)
    checkIdentical(as.matrix(Dr2), ref) 

    # Works with >2 matrices.
    Dr3 <- rbind(Dx, Dy, Dx, Dy)
    ref <- as.matrix(rbind(X, Y, X, Y))
    dimnames(ref) <- NULL
    checkIdentical(as.matrix(Dr3), ref) 
}

test_native_trans <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)

    ref <- as.matrix(t(X))
    dimnames(ref) <- NULL
    checkIdentical(as.matrix(t(Dx)), ref) 

    checkIdentical(as.matrix(t(t(Dx))), as.matrix(X))
}

test_native_ops <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)

    checkIdentical(as.matrix(X+1), as.matrix(Dx+1))
    checkIdentical(as.matrix(X*2), as.matrix(Dx*2))

    ref <- as.matrix(X > 0)
    dimnames(ref) <- list(NULL, NULL)
    checkIdentical(ref, as.matrix(Dx > 0))

    Y <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dy <- DelayedArray(Y)

    checkIdentical(as.matrix(X + Y), as.matrix(Dx + Dy))
    checkIdentical(as.matrix(X * Y), as.matrix(Dx * Dy))
}

test_native_mult <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)

    mat <- matrix(runif(500), ncol=5)
    checkIdentical(unname(as.matrix(X %*% mat)), as.matrix(Dx %*% mat))

    mat2 <- matrix(runif(5000), nrow=5)
    checkIdentical(unname(as.matrix(mat2 %*% X)), as.matrix(mat2 %*% Dx))

    checkIdentical(unname(as.matrix(crossprod(X, t(mat2)))), as.matrix(crossprod(Dx, t(mat2))))
    checkIdentical(unname(as.matrix(crossprod(t(mat2), X))), as.matrix(crossprod(t(mat2), Dx)))

    checkIdentical(unname(as.matrix(tcrossprod(X, t(mat)))), as.matrix(tcrossprod(Dx, t(mat))))
    checkIdentical(unname(as.matrix(tcrossprod(t(mat), X))), as.matrix(tcrossprod(t(mat), Dx)))

    checkIdentical(unname(as.matrix(crossprod(X))), as.matrix(crossprod(Dx)))
    checkIdentical(unname(as.matrix(tcrossprod(X))), as.matrix(tcrossprod(Dx)))
}

test_native_stats <- function() {
    X <- Matrix::rsparsematrix(1000, 100, 0.1)
    Dx <- DelayedArray(X)

    checkEquals(rowSums(X), rowSums(Dx))
    checkEquals(rowMeans(X), rowMeans(Dx))
    checkEquals(colSums(X), colSums(Dx))
    checkEquals(colMeans(X), colMeans(Dx))

    checkEquals(rowMaxs(as.matrix(X)), rowMaxs(Dx))
    checkEquals(rowMins(as.matrix(X)), rowMins(Dx))
    checkEquals(rowRanges(as.matrix(X)), rowRanges(Dx))

    checkEquals(colMaxs(as.matrix(X)), colMaxs(Dx))
    checkEquals(colMins(as.matrix(X)), colMins(Dx))
    checkEquals(colRanges(as.matrix(X)), colRanges(Dx))
}
