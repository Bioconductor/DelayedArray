
test_RleArray <- function()
{
    data <- Rle(1:200000, 125)
    A1 <- RleArray(data, c(62500, 400))
    A2 <- RleArray(data, c(62500, 400), chunksize=1e8)

    on.exit(suppressMessages(setAutoBlockSize()))
    suppressMessages(setAutoBlockSize(10e6))
    rs1 <- rowSums(A1)
    rs2 <- rowSums(A2)
    checkIdentical(rs1, rs2)
    cs1 <- colSums(A1)
    cs2 <- colSums(A2)
    checkIdentical(cs1, cs2)

    ## TODO: Add more tests...
}

test_long_RleArray <- function()
{
    data <- list(Rle(1:5000, 200000), Rle(1:20000, 12500), Rle(1:50000, 20000))
    A <- RleArray(data, dim=c(30000, 75000))
    checkTrue(validObject(A, complete=TRUE))
    checkTrue(is(seed(A), "ChunkedRleArraySeed"))

    ## TODO: Add more tests...
}

test_coercion_to_RleArray <- function()
{
    TYPES <- c("logical", "integer", "double", "complex", "character", "raw")
    for (type in TYPES) {
        x0 <- vector(type)
        from <- as(Rle(x0), "CompressedRleList")
        current <- as(from, "RleArray")
        checkTrue(is(current, "RleMatrix"))
        checkTrue(validObject(current, complete=TRUE))
        checkIdentical(c(0L, 0L), dim(current))
        checkIdentical(type, type(current))
        checkIdentical(matrix(x0, nrow=0, ncol=0), as.matrix(current))
    }

    from <- RleList(compress=FALSE)
    current <- as(from, "RleArray")
    checkTrue(is(current, "RleMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(c(0L, 0L), dim(current))
    checkIdentical("logical", type(current))
    checkIdentical(matrix(nrow=0, ncol=0), as.matrix(current))

    from <- RleList(A=Rle(1, 10), B=Rle(2, 10))
    current <- as(from, "RleArray")
    checkTrue(is(current, "RleMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(c(10L, 2L), dim(current))
    checkIdentical("double", type(current))
    checkIdentical(cbind(A=rep(1, 10), B=rep(2, 10)), as.matrix(current))

    from <- RleList(A=Rle(1, 10), B=Rle(2, 9))
    checkException(as(from, "RleArray"))
}

