test_ConstantArray <- function()
{
    A1 <- ConstantArray(c(500, 200), value=1)
    checkIdentical(as.matrix(A1), matrix(1, 500, 200))
    checkIdentical(as.matrix(A1[1:10,]), matrix(1, 10, 200))
    checkIdentical(as.matrix(A1[,1:10]), matrix(1, 500, 10))
    checkIdentical(as.character(class(A1)), "ConstantMatrix")
    checkTrue(!is_sparse(A1))

    A2 <- ConstantArray(c(500, 20, 10), value=NA)
    checkIdentical(as.array(A2), array(NA, c(500, 20, 10)))
    checkIdentical(as.character(class(A2)), "ConstantArray")
    checkTrue(!is_sparse(A2))

    A3 <- ConstantArray(c(100, 200), value=0)
    checkTrue(is_sparse(A3))
    checkIdentical(as.matrix(A3), matrix(0, 100, 200))

    out <- extract_sparse_array(A3, list(1:10, 1:20))
    checkIdentical(dim(out), c(10L, 20L))
    checkIdentical(nzdata(out), numeric(0))
}
