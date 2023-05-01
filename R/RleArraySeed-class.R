### =========================================================================
### RleArraySeed objects
### -------------------------------------------------------------------------


setClass("RleArraySeed",
    contains="Array",
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
### implements realization of an array-like object as an RleArray object.
### Once writing the array data to the RleRealizationSink object is complete,
### the object will be turned into a ChunkedRleArraySeed object that will be
### used as the seed of the RleArray object.
setClass("RleRealizationSink",
    contains=c("RleArraySeed", "RealizationSink"),
    representation(
        type="character",                     # Single string.
        ## TODO: Add the 2 slots below to make RleRealizationSink
        ## support a RegularArrayGrid of chunks.
        #chunk_grid="RegularArrayGrid",        # Of length N.
        #chunk_runs_along_last_dim="logical",  # Of length N.
        chunks="environment"                  # Of length N (once all the
                                              # chunks are written).
    )
)

setMethod("type", "RleRealizationSink", function(x) x@type)

#setMethod("chunkdim", "RleRealizationSink",
#    function(x) dim(x@chunk_grid[[1L]])
#)

### We support long ChunkedRleArraySeed objects but the chunks cannot be long.
### Note that supporting long chunks would require (at least) that:
###   1) we support long ArrayViewport objects,
###   2) S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a numeric
###      vector.
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

## TODO: Replace ChunkedRleArraySeed above definition with definition below
## to make ChunkedRleArraySeed support a RegularArrayGrid of chunks.
#setClass("ChunkedRleArraySeed", contains="RleRealizationSink")

#setMethod("vertical_slot_names", "ChunkedRleArraySeed",
#    function(x) c("chunk_grid", "chunk_runs_along_last_dim", "chunks")
#)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level chunk accessors
###

.get_chunk <- function(envir, k)
{
    name <- sprintf("%06d", k)
    stopifnot(nchar(name) == 6L)
    get(name, envir=envir, inherits=FALSE)
}

get_chunk_lengths <- function(envir)
{
    ## Too bad we can't just do 'lengths(envir)' for this.
    ## Also would have been nice to be able to just do
    ## 'unlist(eapply(envir, length))' but the list returned by eapply()
    ## is not guaranteed to be sorted and eapply() does not have a 'sorted'
    ## argument. So would need to manually sort it.
    ## Another possibility would be to vapply() on the sorted symbols
    ## returned by 'ls(envir, sorted=TRUE)'.
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
### Validity
###

.validate_RleArraySeed <- function(x)
{
    msg <- S4Arrays:::validate_dim_slot(x, "DIM")
    if (!isTRUE(msg))
        return(msg)
    msg <- S4Arrays:::validate_dimnames_slot(x, x@DIM, "DIMNAMES")
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("RleArraySeed", .validate_RleArraySeed)

.validate_SolidRleArraySeed <- function(x)
{
    ## 'rle' slot.
    if (!is(x@rle, "Rle"))
        return("'rle' slot must be an Rle object")
    x_len <- length(x)
    data_len <- length(x@rle)
    if (x_len != data_len)
        return(paste0("object dimensions [product ", x_len, "] do not ",
                      "match the length of its data [" , data_len, "]"))
    ## Until S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a
    ## numeric vector, we cannot support long SolidRleArraySeed objects.
    if (x_len > .Machine$integer.max)
        return("long SolidRleArraySeed objects are not supported yet")
    TRUE
}
setValidity2("SolidRleArraySeed", .validate_SolidRleArraySeed)

.validate_RleRealizationSink <- function(x)
{
    ## 'type' slot.
    if (!isSingleString(x@type))
        return("'type' slot must be a single string")

    ## TODO: Add the 2 checks below when RleRealizationSink supports a
    ## RegularArrayGrid of chunks.
    ## 'chunk_grid' slot.
    #if (!is(x@chunk_grid, "RegularArrayGrid"))
    #    return("'chunk_grid' slot must be a RegularArrayGrid object")
    #if (!identical(x@DIM, refdim(x@chunk_grid)))
    #    return("'chunk_grid' slot must be a grid that fits the object")
    ## 'chunk_runs_along_last_dim' slot.
    #if (!is.logical(x@chunk_runs_along_last_dim))
    #    return("'chunk_runs_along_last_dim' slot must be a logical vector")

    ## 'chunks' slot.
    if (!is.environment(x@chunks))
        return("'chunks' slot must be an environment")
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
        return(paste0("'x@breakpoints' must be a numeric vector containing ",
                      "strictly sorted positive values"))
    x_len <- length(x)
    data_len <- .get_data_length_from_breakpoints(x@breakpoints)
    if (data_len != x_len)
        return(paste0("length of object data [" , data_len, "] does not ",
                      "match object dimensions [product ", x_len, "]"))
    chunk_lens <- diff(c(0, x@breakpoints))  # chunk lengths as inferred
                                             # from 'breakpoints'
    ## Until S4Vectors:::extract_positions_from_Rle() accepts 'pos' as a
    ## numeric vector, the chunks cannot be long Rle objects.
    if (any(chunk_lens > .Machine$integer.max))
        return("ChunkedRleArraySeed objects do not support long chunks yet")
    # TODO: Check that the chunk lengths as inferred from 'breakpoints'
    # actually match the real ones.
    TRUE
}

### TODO: Replace validity method above with simpler method below when
### ChunkedRleArraySeed supports a RegularArrayGrid of chunks.
#.validate_ChunkedRleArraySeed <- function(x)
#{
#    ## 'chunk_runs_along_last_dim' slot.
#    if (anyNA(x@chunk_runs_along_last_dim))
#        return(paste0("'chunk_runs_along_last_dim' slot must ",
#                      "be a logical vector with no NAs"))
#    ## 'chunks' slot.
#    if (!identical(lengths(x@chunk_grid), get_chunk_lengths(x@chunks)))
#        return("chunk lengths don't match chunking grid element lengths")
#    TRUE
#}
setValidity2("ChunkedRleArraySeed", .validate_ChunkedRleArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("dim", "RleArraySeed", function(x) x@DIM)

setMethod("dimnames", "RleArraySeed",
    function(x) S4Arrays:::simplify_NULL_dimnames(x@DIMNAMES)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from RleArraySeed to Rle or RleList
###

setAs("SolidRleArraySeed", "Rle", function(from) from@rle)

### Return an unnamed RleList object with 1 list element per chunk in 'from'.
setAs("RleRealizationSink", "RleList",
    function(from)
    {
        ## When the object is empty, return a CompressedRleList, NOT a
        ## SimpleRleList object, so unlist() can be used on it.
        if (length(from@chunks) == 0L)
            return(as(Rle(vector(from@type)), "CompressedRleList"))
        chunks <- unname(as.list(from@chunks, sorted=TRUE))
        chunks <- lapply(chunks,
               function(chunk) {
                   chunk_vals <- runValue(chunk)
                   if (typeof(chunk_vals) != from@type) {
                       storage.mode(chunk_vals) <- from@type
                       runValue(chunk) <- chunk_vals
                   }
                   chunk
               })
        as(chunks, "SimpleRleList")
    }
)

### In practice this coercion is not used on an RleRealizationSink instance
### but on a *ChunkedRleArraySeed* instance (e.g. by coercion methods from
### ChunkedRleArraySeed to SolidRleArraySeed or from RleArray to Rle defined
### below in this file).
setAs("RleRealizationSink", "Rle",
    function(from) unlist(as(from, "RleList"), use.names=FALSE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_SolidRleArraySeed <- function(x, index)
{
    x_dim <- dim(x)
    i <- S4Arrays:::to_linear_index(index, x_dim)
    ans <- S4Vectors:::extract_positions_from_Rle(x@rle, i, decoded=TRUE)
    S4Arrays:::set_dim(ans, S4Arrays:::get_Nindex_lengths(index, x_dim))
}
setMethod("extract_array", "SolidRleArraySeed",
    .extract_array_from_SolidRleArraySeed
)

.extract_array_from_ChunkedRleArraySeed <- function(x, index)
{
    x_dim <- dim(x)
    i <- S4Arrays:::to_linear_index(index, x_dim)
    ans <- vector(x@type)
    if (length(i) != 0L) {
        part_idx <- S4Arrays:::get_part_index(i, x@breakpoints)
        split_part_idx <- S4Arrays:::split_part_index(part_idx,
                                                      length(x@breakpoints))
        chunk_idx <- which(lengths(split_part_idx) != 0L)  # chunks to visit
        res <- lapply(chunk_idx, function(i1) {
            chunk <- .get_chunk(x@chunks, i1)
            ## Because a valid ChunkedRleArraySeed object is guaranteed to not
            ## contain long chunks at the moment, 'i2' can be represented as
            ## an integer vector.
            i2 <- as.integer(split_part_idx[[i1]])
            S4Vectors:::extract_positions_from_Rle(chunk, i2, decoded=TRUE)
        })
        res <- c(list(ans), res)
        ans <- unlist(res, use.names=FALSE)[S4Arrays:::get_rev_index(part_idx)]
    }
    S4Arrays:::set_dim(ans, S4Arrays:::get_Nindex_lengths(index, x_dim))
}
setMethod("extract_array", "ChunkedRleArraySeed",
    .extract_array_from_ChunkedRleArraySeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleRealizationSink constructor
###

### NOT exported!
RleRealizationSink <- function(dim, dimnames=NULL, type="double")
{
    dim <- S4Arrays:::normarg_dim(dim)
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    chunks <- new.env(hash=TRUE, parent=emptyenv())
    new2("RleRealizationSink", DIM=dim, DIMNAMES=dimnames,
                               type=type, chunks=chunks)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SolidRleArraySeed and ChunkedRleArraySeed low-level constructors
###

.new_SolidRleArraySeed <- function(data, dim, dimnames)
{
    if (!is(data, "Rle"))
        stop(wmsg("invalid 'data'"))
    new2("SolidRleArraySeed", DIM=dim, DIMNAMES=dimnames, rle=data)
}

append_Rle_to_sink <- function(x, sink)
{
    stopifnot(is(x, "Rle"))
    if (length(x) == 0L)
        return()  # nothing to do
    if (sink@type == "integer") {
        x_vals <- runValue(x)
        ## Replace integer-Rle with raw-Rle if this doesn't loose
        ## information.
        if (!S4Vectors:::anyMissingOrOutside(x_vals, 0L, 255L))
            runValue(x) <- as.raw(x_vals)
    }
    .append_chunk(sink@chunks, x)
}

### This coercion is used by .new_ChunkedRleArraySeed() and by the
### coercion method from RleRealizationSink to RleArray.
setAs("RleRealizationSink", "ChunkedRleArraySeed",
    function(from)
    {
        breakpoints <- cumsum(as.double(get_chunk_lengths(from@chunks)))
        new2("ChunkedRleArraySeed", from, breakpoints=breakpoints)
    }
)

### 'data' is assumed to be a list-like object where all the list elements
### are Rle objects of type 'type'. This is NOT checked!
.new_ChunkedRleArraySeed <- function(data, dim, dimnames, type)
{
    sink <- RleRealizationSink(dim, dimnames, type)
    for (k in seq_along(data))
        append_Rle_to_sink(data[[k]], sink)
    as(sink, "ChunkedRleArraySeed")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleArraySeed high-level constructor
###

### 'data' must be an Rle object. This is not checked.
### If 'chunksize' is not specified, return a SolidRleArraySeed instance.
### Otherwise return a ChunkedRleArraySeed instance.
.make_RleArraySeed_from_Rle <- function(data, dim, dimnames, chunksize=NULL)
{
    dim <- S4Arrays:::normarg_dim(dim)
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    ans_len <- length(data)
    if (ans_len != prod(dim))
        stop(wmsg("length of input data [" , ans_len, "] does not ",
                  "match object dimensions [product ", prod(dim), "]"))
    if (ans_len > .Machine$integer.max)
        stop(wmsg("Input data is too long (>= 2^31). Please supply an ",
                  "ordinary list of Rle objects instead, or an RleList ",
                  "object, or a DataFrame object where all the columns ",
                  "are Rle objects."))
    if (is.null(chunksize))
        return(.new_SolidRleArraySeed(data, dim, dimnames))
    ans_type <- typeof(runValue(data))
    ## FIXME: breakInChunks() does not accept a 'totalsize' >= 2^31 at
    ## the moment so this won't work on a long Rle.
    partitioning <- breakInChunks(ans_len, chunksize=chunksize)
    data <- relist(data, partitioning)
    .new_ChunkedRleArraySeed(data, dim, dimnames, ans_type)
}

### 'data' is assumed to be a list-like object.
.infer_dim <- function(data)
{
    data_dim <- dim(data)
    if (!is.null(data_dim)) {
        ## Potecting against objects with a weird dim().
        if (!is.integer(data_dim))
            stop(wmsg("Please specify the 'dim' argument. ",
                      "'dim(data)' is not an integer vector so ",
                      "cannot be used as a fallback for 'dim'."))
        return(data_dim)
    }
    data_ncol <- length(data)
    if (data_ncol == 0L)
        return(c(0L, 0L))
    data_lens <- lengths(data, use.names=FALSE)
    data_nrow <- data_lens[[1L]]
    if (!all(data_lens == data_nrow))
        stop(wmsg("Please specify the 'dim' argument. ",
                  "The dimensions of the object to ",
                  "construct cannot be inferred from 'data'."))
    c(data_nrow, data_ncol)
}

.infer_dimnames <- function(data)
{
    data_dimnames <- dimnames(data)
    if (!is.null(data_dimnames))
        return(data_dimnames)
    data_names <- names(data)
    if (!is.null(data_names))
        return(list(NULL, data_names))
    NULL
}

### 'data' must be a list-like object where all the list elements are
### Rle objects.
.normalize_data_as_list_of_Rles <- function(data)
{
    if (length(data) == 0L)
        return(list())
    if (is(data, "CompressedRleList"))
        return(as.list(data))
    if (!is(data, "RleList")) {
        ok <- vapply(data, is, logical(1), "Rle", USE.NAMES=FALSE)
        if (!all(ok))
            stop(wmsg("all the list elements in the input object ",
                      "must be Rle objects"))
    }
    ## Turn 'data' into an ordinary list where all the list elements are
    ## Rle objects of the same type.
    data0 <- lapply(data, function(x) runValue(x)[integer(0)])
    data_type <- typeof(unlist(data0, use.names=FALSE))
    lapply(data,
        function(x) {
            x_vals <- runValue(x)
            if (typeof(x_vals) != data_type) {
                storage.mode(x_vals) <- data_type
                runValue(x) <- x_vals
            }
            x
        })
}

### 'data' must be a list-like object where all the list elements are Rle
### objects. If 'data' has length 0, return a SolidRleArraySeed instance.
### Otherwise return a ChunkedRleArraySeed instance.
.make_RleArraySeed_from_list_or_List <- function(data, dim, dimnames)
{
    if (missing(dim)) {
        dim <- .infer_dim(data)
        if (missing(dimnames))
            dimnames <- .infer_dimnames(data)
    } else {
        dim <- S4Arrays:::normarg_dim(dim)
        ans_len <- sum(lengths(data))  # can be >= 2^31
        if (ans_len != prod(dim))
            stop(wmsg("total length of input data [" , ans_len, "] does not ",
                      "match object dimensions [product ", prod(dim), "]"))
    }
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    if (length(data) == 0L) {
        unlisted_data <- unlist(data, use.names=FALSE)
        if (is.null(unlisted_data))
            unlisted_data <- Rle()
        return(.new_SolidRleArraySeed(unlisted_data, dim, dimnames))
    }
    data <- .normalize_data_as_list_of_Rles(data)
    ans_type <- typeof(runValue(data[[1L]]))
    .new_ChunkedRleArraySeed(data, dim, dimnames, ans_type)
}

### 'data' must be an Rle object or a list-like object where all the list
### elements are Rle objects. In the former case, and if 'chunksize' is
### not specified, a SolidRleArraySeed object is returned. Otherwise a
### ChunkedRleArraySeed object is returned.
### NOT exported!
RleArraySeed <- function(data, dim, dimnames, chunksize=NULL)
{
    if (is(data, "Rle"))
        return(.make_RleArraySeed_from_Rle(data, dim, dimnames, chunksize))
    if (is(data, "list_OR_List")) {
        if (!is.null(chunksize))
            warning(wmsg("'chunksize' is currently ignored ",
                         "when the input is a list-like object"))
        return(.make_RleArraySeed_from_list_or_List(data, dim, dimnames))
    }
    stop(wmsg("'data' must be an Rle object or a list-like object ",
              "where all the list elements are Rle objects"))
}

setAs("ChunkedRleArraySeed", "SolidRleArraySeed",
    function(from) RleArraySeed(as(from, "Rle"), dim(from), dimnames(from))
)

