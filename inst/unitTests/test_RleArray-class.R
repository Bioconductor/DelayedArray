DEFAULT_BLOCK_SIZE <- DelayedArray:::DEFAULT_BLOCK_SIZE

test_RleArray <- function()
{
    rle <- Rle(1:200000, 125)
    A1 <- RleArray(rle, c(62500, 400))
    A2 <- RleArray(rle, c(62500, 400), chunksize=1e8)

    on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
    options(DelayedArray.block.size=10e6)
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
    ## Right now it's not possible to create a long RleArray object with
    ## the RleArray() constructor function. So we use the low-level RleArray
    ## construction API to do this:
    RleRealizationSink <- DelayedArray:::RleRealizationSink
    append_Rle_to_sink <- DelayedArray:::.append_Rle_to_sink
    sink <- RleRealizationSink(c(30000L, 75000L), type="integer")
    #rle1 <- Rle(1:500000, 2000)
    rle1 <- Rle(1:5000, 200000)
    append_Rle_to_sink(rle1, sink)
    #rle2 <- Rle(1:2000000, 125)
    rle2 <- Rle(1:20000, 12500)
    append_Rle_to_sink(rle2, sink)
    #rle3 <- Rle(1:5000000, 200)
    rle3 <- Rle(1:50000, 20000)
    append_Rle_to_sink(rle3, sink)
    A <- as(sink, "RleArray")

    checkTrue(validObject(A, complete=TRUE))
    checkTrue(is(seed(A), "ChunkedRleArraySeed"))

    ## TODO: Add more tests...
}

