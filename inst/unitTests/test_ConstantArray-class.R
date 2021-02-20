test_ConstantArray <- function()
{
    A1 <- ConstantArray(c(500, 200), value=1)
    checkIdentical(as.matrix(A1), matrix(1, 500, 200))
    checkIdentical(as.matrix(A1[1:10,]), matrix(1, 10, 200))
    checkIdentical(as.matrix(A1[,1:10]), matrix(1, 500, 10))

    checkIdentical(as.character(class(A1)), "ConstantMatrix")
    checkTrue(is(A1, "ConstantArray"))
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

test_ConstantArray_other <- function()
{
    # Testing some of the more odd types we can put in here.
    Ac <- ConstantArray(c(500, 200), value="Aaron")
    checkIdentical(as.matrix(Ac), matrix("Aaron", 500, 200))
    checkTrue(!is_sparse(Ac))

    Acs <- ConstantArray(c(500, 200), value="")
    checkIdentical(as.matrix(Acs), matrix("", 500, 200))
    checkTrue(is_sparse(Acs))

    Al <- ConstantArray(c(500, 200), value=list("Aaron"))
    checkIdentical(as.matrix(Al), matrix(list("Aaron"), 500, 200))

    Al2 <- ConstantArray(c(500, 200), value=list(letters))
    checkIdentical(as.matrix(Al2), matrix(list(letters), 500, 200))

    # Checking error states.
    checkException(ConstantArray(c(500, -1), value=0), silent=TRUE)
    checkException(ConstantArray(c(500, 200), value=letters), silent=TRUE)
}

test_ConstantArray_coercion <- function()
{
    A1 <- ConstantArray(c(500, 200), value=1)
    checkIdentical(A1, as(A1, "ConstantArray"))

    seed <- ConstantArraySeed(c(500, 200), value=1)
    A2 <- new("ConstantArray", seed=seed)
    checkIdentical(as.character(class(A2)[1]), "ConstantArray")
    A2m <- as(A2, "ConstantMatrix")
    checkIdentical(A2m, DelayedArray(seed))
}
