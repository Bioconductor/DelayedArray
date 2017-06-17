### =========================================================================
### RleArray objects
### -------------------------------------------------------------------------


setClass("RleArraySeed",
    representation(
        "VIRTUAL",
        ## Must use upper case or won't be able to extend the class.
        ## See https://stat.ethz.ch/pipermail/r-devel/2017-June/074383.html
        DIM="integer",
        DIMNAMES="list"
    )
)

### We don't support long SolidRleArraySeed objects yet! This would first
### require that S4Vectors:::extract_positions_from_Rle() accepts 'pos' as
### a numeric vector.
setClass("SolidRleArraySeed",
    contains="RleArraySeed",
    representation(
        rle="Rle"
    )
)

### The RleRealizationSink class is a concrete RealizationSink subclass that
### implements realization of an array-like object as an RleArray object that
### will have a ChunkedRleArraySeed seed (once writting to the sink is
### complete).
setClass("RleRealizationSink",
    contains=c("RleArraySeed", "RealizationSink"),
    representation(
        type="character",
        chunks="environment"
    )
)

### We support long ChunkedRleArraySeed objects but for now the chunks
### cannot be long. Supporting long chunks would require that
### S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a numeric
### vector.
setClass("ChunkedRleArraySeed",
    contains="RleRealizationSink",
    representation(
        ## A numeric vector of length the nb of chunks. Contains the cumulated
        ## lengths of the chunks so must be "numeric" (and not "integer") to
        ## support long objects. A chunk cannot be empty so 'breakpoints' must
        ## contain *strictly* sorted positive values.
        ## If the object is of length 0, then 'breakpoints' is empty.
        ## Otherwise, its last element must equal the length of the object.
        breakpoints="numeric"
    )
)

.validate_RleArraySeed <- function(x)
{
    msg <- validate_dim_slot(x, "DIM")
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dimnames_slot(x, x@DIM, "DIMNAMES")
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("RleArraySeed", .validate_RleArraySeed)

.validate_SolidRleArraySeed <- function(x)
{
    ## 'rle' slot.
    if (!is(x@rle, "Rle"))
        return(wmsg2("'x@rle' must be an Rle object"))
    x_len <- length(x)
    data_len <- length(x@rle)
    if (x_len != data_len)
        return(wmsg2("object dimensions [product ", x_len, "] do not ",
                     "match the length of its data [" , data_len, "]"))
    ## Until S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a
    ## numeric vector, we cannot support long SolidRleArraySeed objects.
    if (x_len > .Machine$integer.max)
        return(wmsg2("long SolidRleArraySeed objects are not supported yet"))
    TRUE
}
setValidity2("SolidRleArraySeed", .validate_SolidRleArraySeed)

.validate_RleRealizationSink <- function(x)
{
    ## 'type' slot.
    if (!isSingleString(x@type))
        return(wmsg2("'x@type' must be a single string"))
    ## 'chunks' slot.
    if (!is.environment(x@chunks))
        return(wmsg2("'x@chunks' must be an environment"))
    # TODO: Validate the content of 'chunks'.
    TRUE
}
setValidity2("RleRealizationSink", .validate_RleRealizationSink)

.get_data_length_from_breakpoints <- function(breakpoints)
{
    breakpoints_len <- length(breakpoints)
    if (breakpoints_len == 0L) 0L else breakpoints[[breakpoints_len]]
}

.validate_ChunkedRleArraySeed <- function(x)
{
    ## 'breakpoints' slot.
    if (!is.numeric(x@breakpoints)
     || S4Vectors:::anyMissing(x@breakpoints)
     || is.unsorted(x@breakpoints, strictly=TRUE)
     || length(x@breakpoints) != 0L && x@breakpoints[[1L]] <= 0L)
        return(wmsg2("'x@breakpoints' must be a numeric vector containing ",
                     "strictly sorted positive values"))
    x_len <- length(x)
    data_len <- .get_data_length_from_breakpoints(x@breakpoints)
    if (data_len != x_len)
        return(wmsg2("length of object data [" , data_len, "] does not ",
                     "match object dimensions [product ", x_len, "]"))
    chunk_lens <- diff(c(0, x@breakpoints))  # chunk lengths as inferred from
                                             # 'breakpoints'
    ## Until S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a
    ## numeric vector, the chunks cannot be long Rle objects.
    if (any(chunk_lens > .Machine$integer.max))
        return(wmsg2("ChunkedRleArraySeed objects do not support ",
                     "long chunks yet"))
    # TODO: Check that the chunk lengths as inferred from 'breakpoints'
    # actually match the real ones.
    TRUE
}
setValidity2("ChunkedRleArraySeed", .validate_ChunkedRleArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("dim", "RleArraySeed", function(x) x@DIM)

setMethod("length", "RleArraySeed", function(x) prod(dim(x)))

setMethod("dimnames", "RleArraySeed",
    function(x)
    {
        ans <- x@DIMNAMES
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to Rle objects
###

setAs("SolidRleArraySeed", "Rle", function(from) from@rle)

### In practice this coercion is used on a *ChunkedRleArraySeed* object
### by the coercion method from ChunkedRleArraySeed to SolidRleArraySeed
### defined in this file below.
setAs("RleRealizationSink", "Rle",
    function(from)
    {
        ans <- Rle(match.fun(from@type)(0))
        if (length(from@chunks) == 0L)
            return(ans)
        list_of_Rles <- c(list(ans), unname(as.list(from@chunks, sorted=TRUE)))
        do.call("c", list_of_Rles)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level chunk accessors
###

.get_chunk <- function(envir, k)
{
    name <- sprintf("%06d", k)
    stopifnot(nchar(name) == 6L)
    get(name, envir=envir, inherits=FALSE)
}

.get_chunk_lens <- function(envir)
{
    ## Too bad we can't just do 'lengths(envir)' for this.
    ## Also would have been nice to be able to just do
    ## 'unlist(eapply(envir, length))' but the list returned by eapply()
    ## is not guaranteed to be sorted and eapply() does not have a 'sorted'
    ## argument. So would need to manually sort it.
    ## Another possibility would be to vapply() on the sorted symbols returned
    ## by 'ls(envir, sorted=TRUE)'.
    vapply(seq_len(length(envir)),
           function(k) length(.get_chunk(envir, k)),
           numeric(1))
}

.set_chunk <- function(envir, k, chunk)
{
    name <- sprintf("%06d", k)
    stopifnot(nchar(name) == 6L)
    assign(name, chunk, envir=envir)
}

.append_chunk <- function(envir, chunk)
{
    .set_chunk(envir, length(envir) + 1L, chunk)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()
###

.subset_SolidRleArraySeed_as_array <- function(seed, index)
{
    seed_dim <- dim(seed)
    i <- to_linear_index(index, seed_dim)
    ans <- S4Vectors:::extract_positions_from_Rle(seed@rle, i, decoded=TRUE)
    dim(ans) <- get_Nindex_lengths(index, seed_dim)
    ans
}
setMethod("subset_seed_as_array", "SolidRleArraySeed",
    .subset_SolidRleArraySeed_as_array
)

.subset_ChunkedRleArraySeed_as_array <- function(seed, index)
{
    seed_dim <- dim(seed)
    i <- to_linear_index(index, seed_dim)
    ans <- match.fun(seed@type)(0)
    if (length(i) != 0L) {
        part_idx <- get_part_index(i, seed@breakpoints)
        split_part_idx <- split_part_index(part_idx, length(seed@breakpoints))
        chunk_idx <- which(lengths(split_part_idx) != 0L)  # chunks to visit
        res <- lapply(chunk_idx, function(i1) {
            chunk <- .get_chunk(seed@chunks, i1)
            ## Because a valid ChunkedRleArraySeed object is guaranteed to not
            ## contain long chunks at the moment, 'i2' can be represented as
            ## an integer vector.
            i2 <- as.integer(split_part_idx[[i1]])
            S4Vectors:::extract_positions_from_Rle(chunk, i2, decoded=TRUE)
        })
        res <- c(list(ans), res)
        ans <- unlist(res, use.names=FALSE)[get_rev_index(part_idx)]
    }
    dim(ans) <- get_Nindex_lengths(index, seed_dim)
    ans
}
setMethod("subset_seed_as_array", "ChunkedRleArraySeed",
    .subset_ChunkedRleArraySeed_as_array
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Construction of RleRealizationSink and RleArraySeed objects
###

### NOT exported!
RleRealizationSink <- function(dim, dimnames=NULL, type="double")
{
    if (is.null(dimnames))
        dimnames <- vector("list", length(dim))
    chunks <- new.env(hash=TRUE, parent=emptyenv())
    new2("RleRealizationSink", DIM=dim, DIMNAMES=dimnames,
                               type=type, chunks=chunks)
}

### Ignores 'viewport'.
setMethod("write_to_sink", c("Rle", "RleRealizationSink"),
    function(x, sink, viewport)
    {
        if (length(x) == 0L)
            return()  # nothing to do
        if (sink@type == "integer") {
            run_values <- runValue(x)
            ## Replace integer-Rle with raw-Rle if this doesn't loose
            ## information.
            if (!S4Vectors:::anyMissingOrOutside(run_values, 0L, 255L))
                runValue(x) <- as.raw(run_values)
        }
        .append_chunk(sink@chunks, x)
    }
)

### This coercion is used by the RleArraySeed() constructor and by the
### coercion method from RleRealizationSink to RleArray.
setAs("RleRealizationSink", "ChunkedRleArraySeed",
    function(from)
    {
        breakpoints <- cumsum(as.double(.get_chunk_lens(from@chunks)))
        new2("ChunkedRleArraySeed", from, breakpoints=breakpoints)
    }
)

### NOT exported!
RleArraySeed <- function(rle, dim, dimnames=NULL, chunksize=NULL)
{
    if (!is.numeric(dim))
        stop(wmsg("the supplied dim vector must be numeric"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(dimnames))
        dimnames <- vector("list", length=length(dim))

    if (is.null(chunksize))
        return(new2("SolidRleArraySeed", DIM=dim, DIMNAMES=dimnames, rle=rle))

    type <- typeof(runValue(rle))
    sink <- RleRealizationSink(dim, dimnames, type)
    ## FIXME: breakInChunks() does not accept a 'totalsize' >= 2^31 at the
    ## moment so this won't work on a long Rle.
    partitioning <- breakInChunks(length(rle), chunksize)
    rle_list <- relist(rle, partitioning)
    for (k in seq_along(rle_list))
        write_to_sink(rle_list[[k]], sink)
    as(sink, "ChunkedRleArraySeed")
}

setAs("ChunkedRleArraySeed", "SolidRleArraySeed",
    function(from) RleArraySeed(as(from, "Rle"), dim(from), dimnames(from))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleArray and RleMatrix objects
###

setClass("RleArray", contains="DelayedArray")

setClass("RleMatrix", contains=c("DelayedMatrix", "RleArray"))

### Automatic coercion method from RleArray to RleMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
setAs("RleArray", "RleMatrix", function(from) new2("RleMatrix", from))

### For internal use only.
setMethod("matrixClass", "RleArray", function(x) "RleMatrix")

.validate_RleArray <- function(x)
{
    if (!is(seed(x), "RleArraySeed"))
        return(wmsg2("'seed(x)' must be an RleArraySeed object"))
    if (!is_pristine(x))
        return(wmsg2("'x' carries delayed operations on it"))
    TRUE
}

setValidity2("RleArray", .validate_RleArray)

setAs("ANY", "RleMatrix",
    function(from) as(as(from, "RleArray"), "RleMatrix")
)

setMethod("DelayedArray", "RleArraySeed",
    function(seed) new_DelayedArray(seed, Class="RleArray")
)

### Works directly on an RleArraySeed object, in which case it must be called
### with a single argument.
RleArray <- function(rle, dim, dimnames=NULL, chunksize=NULL)
{
    if (is(rle, "RleArraySeed")) {
        if (!(missing(dim) && is.null(dimnames) && is.null(chunksize)))
            stop(wmsg("RleArray() must be called with a single argument ",
                      "when passed an RleArraySeed object"))
        seed <- rle
    } else {
        seed <- RleArraySeed(rle, dim, dimnames=dimnames, chunksize=chunksize)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Realization as an RleArray object
###

setMethod("write_to_sink", c("array", "RleRealizationSink"),
    function(x, sink, viewport)
    {
        if (is.null(viewport)) {
            stopifnot(identical(dim(x), dim(sink)))
        } else {
            stopifnot(identical(dim(x), dim(viewport)),
                      identical(dim(sink), refdim(viewport)))
        }
        ## 'viewport' is ignored!
        write_to_sink(Rle(x), sink)
    }
)

setAs("RleRealizationSink", "RleArray",
    function(from) RleArray(as(from, "ChunkedRleArraySeed"))
)

setAs("RleRealizationSink", "DelayedArray", function(from) as(from, "RleArray"))

.as_RleArray <- function(from)
{
    sink <- RleRealizationSink(dim(from), dimnames(from), type(from))
    write_to_sink(from, sink, NULL)
    as(sink, "RleArray")
}

setAs("ANY", "RleArray", .as_RleArray)

### Automatic coercion methods from DelayedArray to RleArray and from
### DelayedMatrix to RleMatrix silently return broken objects (unfortunately
### these dummy automatic coercion methods don't bother to validate the object
### they return). So we overwrite them.
setAs("DelayedArray", "RleArray", .as_RleArray)
setAs("DelayedMatrix", "RleMatrix", .as_RleArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Switching between DataFrame and RleMatrix representation
###

### From DataFrame to RleMatrix.
.from_DataFrame_to_RleMatrix <- function(from)
{
    as(DelayedArray(from), "RleMatrix")
}

setAs("DataFrame", "RleMatrix", .from_DataFrame_to_RleMatrix)
setAs("DataFrame", "RleArray", .from_DataFrame_to_RleMatrix)

### From RleMatrix to DataFrame.
.from_RleMatrix_to_DataFrame <- function(from)
{
    ## We mangle the colnames exactly like as.data.frame() would do.
    ans_colnames <- colnames(as.data.frame(from[0L, ]))
    rle <- as(seed(from), "Rle")
    partitioning <- PartitioningByEnd(nrow(from) * seq_len(ncol(from)),
                                      names=ans_colnames)
    listData <- as.list(relist(rle, partitioning))
    new2("DataFrame", listData=listData,
                      nrows=nrow(from),
                      rownames=rownames(from))
}

setAs("RleMatrix", "DataFrame", .from_RleMatrix_to_DataFrame)

### From DelayedMatrix to DataFrame.
setAs("DelayedMatrix", "DataFrame",
    function(from) as(as(from, "RleMatrix"), "DataFrame")
)

