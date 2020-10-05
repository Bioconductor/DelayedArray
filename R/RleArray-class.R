### =========================================================================
### RleArray objects
### -------------------------------------------------------------------------
###
### Note that we could just wrap an RleArraySeed object in a DelayedArray
### object to represent and manipulate an Rle-encoded array as a DelayedArray
### object. So, strictly speaking, we don't really need the RleArray and
### RleMatrix classes. However, we define these classes mostly for cosmetic
### reasons, that is, to hide the DelayedArray and DelayedMatrix classes
### from the user. So the user will see and manipulate RleArray and
### RleMatrix objects instead of DelayedArray and DelayedMatrix objects.
###

setClass("RleArray",
    contains="DelayedArray",
    representation(seed="RleArraySeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "RleArraySeed",
    function(seed) new_DelayedArray(seed, Class="RleArray")
)

### 'data' must be an Rle object or a list-like object where all the list
### elements are Rle objects. In the former case, and if 'chunksize' is
### not specified, the returned RleArray will be a "solid" RleArray i.e.
### its seed will be a SolidRleArraySeed object. Otherwise it will be
### "chunked" i.e. its seed will be a ChunkedRleArraySeed object.
### In addition RleArray() also works directly on an RleArraySeed object,
### in which case it must be called with a single argument.
RleArray <- function(data, dim, dimnames, chunksize=NULL)
{
    if (is(data, "RleArraySeed")) {
        if (!(missing(dim) && missing(dimnames) && is.null(chunksize)))
            stop(wmsg("RleArray() must be called with a single argument ",
                      "when passed an RleArraySeed object"))
        seed <- data
    } else {
        seed <- RleArraySeed(data, dim, dimnames, chunksize=chunksize)
    }
    DelayedArray(seed)
}

### Deconstruction.
### The returned Rle has the length of the RleArray object.
setAs("RleArray", "Rle", function(from) as(from@seed, "Rle"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleMatrix objects
###

setClass("RleMatrix", contains=c("RleArray", "DelayedMatrix"))

### Required for DelayedArray internal business.
setMethod("matrixClass", "RleArray", function(x) "RleMatrix")

### Automatic coercion method from RleArray to RleMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
setAs("RleArray", "RleMatrix", function(from) new2("RleMatrix", from))

### The user should not be able to degrade an RleMatrix object to
### an RleArray object so 'as(x, "RleArray", strict=TRUE)' should
### fail or be a no-op when 'x' is an RleMatrix object. Making this
### coercion a no-op seems to be the easiest (and safest) way to go.
setAs("RleMatrix", "RleArray", function(from) from)  # no-op

setAs("ANY", "RleMatrix",
    function(from) as(as(from, "RleArray"), "RleMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Switching between RleList and RleMatrix
###

### From RleList to RleMatrix.
### Return an RleMatrix object with 1 column per list element in the RleList
### object. More precisely, when the RleList object has at least 1 column,
### the returned RleMatrix object is "chunked" and has 1 chunk per column
### (so 1 chunk per list element in the original RleList object).
.from_RleList_to_RleMatrix <- function(from)
{
    if (length(from) != 0L) {
        from_lens <- lengths(from, use.names = FALSE)
        from1_len <- from_lens[[1L]]
        if (!all(from_lens == from1_len))
            stop(wmsg("all the list elements in the RleList object ",
                      "to coerce must have the same length"))
    }
    RleArray(from)
}

setAs("RleList", "RleArray", .from_RleList_to_RleMatrix)

### From RleMatrix to RleList.
### Return an RleList object with 1 list element per column in the RleMatrix
### object. More precisely, when the length of the RleMatrix object is < 2^31,
### then a CompressedRleList object is returned. Otherwise, a SimpleRleList
### object is returned. Note that in this case, the current implementation
### only knows how to handle the situation where the RleMatrix has 1 chunk
### per column.
.from_RleMatrix_to_RleList <- function(from)
{
    if (length(from) <= .Machine$integer.max) {
        unlisted_ans <- as(from, "Rle")
        ans_breakpoints <- seq_len(ncol(from)) * nrow(from)
        ans_partitioning <- PartitioningByEnd(ans_breakpoints,
                                              names=colnames(from))
        return(relist(unlisted_ans, ans_partitioning))
    }
    ## 'from@seed' is guaranteed to be a ChunkedRleArraySeed object.
    chunk_lens <- get_chunk_lengths(from@seed@chunks)
    if (!all(chunk_lens == nrow(from)))
        stop(wmsg("coercing this RleMatrix object to an RleArray object ",
                  "is not supported at the moment"))
    setNames(as(from@seed, "RleList"), colnames(from))
}

setAs("RleMatrix", "RleList", .from_RleMatrix_to_RleList)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Switching between DataFrame and RleMatrix
###

### From DataFrame to RleMatrix.
.from_DataFrame_to_RleMatrix <- function(from)
{
    as(DelayedArray(from), "RleMatrix")
}

setAs("DataFrame", "RleArray", .from_DataFrame_to_RleMatrix)

### From RleMatrix to DataFrame.
.from_RleMatrix_to_DataFrame <- function(from)
{
    listData <- as.list(as(from, "RleList"))
    ans <- S4Vectors:::new_DataFrame(listData, nrows=nrow(from))
    rownames(ans) <- rownames(from)
    ans
}

setAs("RleMatrix", "DataFrame", .from_RleMatrix_to_DataFrame)

### From DelayedMatrix to DataFrame.
setAs("DelayedMatrix", "DataFrame",
    function(from) as(as(from, "RleMatrix"), "DataFrame")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Realization as an RleArray object
###

### WARNING: This method assumes that the blocks are "linear" and being
### written in order. Even though this is still the case, this will change
### soon and this change will break it!
### FIXME: The method below must write the block to 'sink' at the location
### specified by the supplied 'viewport'.
setMethod("write_block", "RleRealizationSink",
    function(sink, viewport, block)
    {
        ## 'viewport' is ignored!
        append_Rle_to_sink(Rle(block), sink)
        sink
    }
)

setAs("RleRealizationSink", "RleArray",
    function(from) DelayedArray(as(from, "ChunkedRleArraySeed"))
)

setAs("RleRealizationSink", "DelayedArray", function(from) as(from, "RleArray"))

.as_RleArray <- function(from)
{
    sink <- RleRealizationSink(dim(from), dimnames(from), type(from))
    BLOCK_write_to_sink(sink, from)
    as(sink, "RleArray")
}

setAs("ANY", "RleArray", .as_RleArray)

### Automatic coercion methods from DelayedArray to RleArray and from
### DelayedMatrix to RleMatrix silently return broken objects (unfortunately
### these dummy automatic coercion methods don't bother to validate the object
### they return). So we overwrite them.
setAs("DelayedArray", "RleArray", .as_RleArray)
setAs("DelayedMatrix", "RleMatrix", .as_RleArray)

