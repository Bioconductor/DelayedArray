### =========================================================================
### RleArray objects
### -------------------------------------------------------------------------


### NOT exported!
setClass("RleArraySeed",
    representation(
        rle="Rle",
        dim="integer",
        dimnames="list"
    )
)

.validate_RleArraySeed <- function(x)
{
    ## 'rle' slot.
    if (!is(x@rle, "Rle"))
        return(wmsg("'x@rle' must be an Rle object"))
    ## 'dim' slot.
    if (!is.integer(x@dim))
        return(wmsg("'x@dim' must be an integer vector"))
    x_ndim <- length(x@dim)
    if (x_ndim == 0L)
        return(wmsg("'x@dim' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@dim, 0L))
        return(wmsg("'x@dim' cannot contain negative or NA values"))
    p <- prod(x@dim)
    if (p != length(x@rle))
        return(wmsg("object dimensions [product ", p, "] do not match ",
                    "its length [" , length(x@rle), "]"))
    ## 'dimnames' slot.
    if (!is.list(x@dimnames))
        return(wmsg("'x@dimnames' must be a list of length"))
    if (length(x@dimnames) != x_ndim)
        return(wmsg("length of 'x@dimnames' must match that of 'x@dim'"))
    if (!all(vapply(seq_len(x_ndim),
                    function(n) {
                      dn <- x@dimnames[[n]]
                      if (is.null(dn))
                          return(TRUE)
                      is.character(dn) && length(dn) == x@dim[[n]]
                    },
                    logical(1),
                    USE.NAMES=FALSE)))
        return(wmsg("every list element in 'x@dimnames' must be NULL or ",
                    "a character vector along the corresponding dimension"))
    TRUE
}

setValidity2("RleArraySeed", .validate_RleArraySeed)

setMethod("dim", "RleArraySeed", function(x) x@dim)

setMethod("length", "RleArraySeed", function(x) prod(dim(x)))

setMethod("dimnames", "RleArraySeed",
    function(x)
    {
        ans <- x@dimnames
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()
###

.subset_RleArraySeed_as_array <- function(seed, index)
{
    seed_dim <- dim(seed)
    i <- to_linear_index(index, seed_dim)
    ans <- decode(seed@rle[i])
    dim(ans) <- get_Nindex_lengths(index, seed_dim)
    ans
}

setMethod("subset_seed_as_array", "RleArraySeed",
    .subset_RleArraySeed_as_array
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleArraySeed constructor
###

### NOT exported!
RleArraySeed <- function(rle, dim, dimnames=NULL)
{
    if (!is.numeric(dim))
        stop(wmsg("the supplied dim vector must be numeric"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(dimnames))
        dimnames <- vector("list", length=length(dim))
    new2("RleArraySeed", rle=rle, dim=dim, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleArray and RleMatrix objects
###

setClass("RleArray", contains="DelayedArray")

setClass("RleMatrix", contains=c("DelayedMatrix", "RleArray"))

### Automatic coercion method from RleArray to RleMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
setAs("RleArray", "RleMatrix", function(from) new("RleMatrix", from))

### For internal use only.
setMethod("matrixClass", "RleArray", function(x) "RleMatrix")

.validate_RleArray <- function(x)
{
    if (!is(seed(x), "RleArraySeed"))
        return(wmsg("'seed(x)' must be an RleArraySeed object"))
    if (!is_pristine(x))
        return(wmsg("'x' carries delayed operations on it"))
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
RleArray <- function(rle, dim, dimnames=NULL)
{
    if (is(rle, "RleArraySeed")) {
        if (!(missing(dim) && is.null(dimnames)))
            stop(wmsg("RleArray() must be called with a single argument ",
                      "when passed an RleArraySeed object"))
        seed <- rle
    } else {
        seed <- RleArraySeed(rle, dim, dimnames=dimnames)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RleRealizationSink objects
###
### The RleRealizationSink class is a concrete RealizationSink subclass that
### implements realization of an array-like object as an RleArray object.
###

setClass("RleRealizationSink",
    contains="RealizationSink",
    representation(
        dim="integer",
        dimnames="list",
        type="character",
        dump="environment"
    )
)

RleRealizationSink <- function(dim, dimnames=NULL, type="double")
{
    if (is.null(dimnames))
        dimnames <- vector("list", length(dim))
    dump <- new.env(parent=emptyenv())
    new("RleRealizationSink", dim=dim, dimnames=dimnames, type=type, dump=dump)
}

setMethod("write_to_sink", c("array", "RleRealizationSink"),
    function(x, sink, offsets=NULL)
    {
        x_dim <- dim(x)
        sink_dim <- sink@dim
        stopifnot(length(x_dim) == length(sink_dim))
        if (is.null(offsets)) {
            stopifnot(identical(x_dim, sink_dim))
        } else {
            stopifnot(length(x_dim) == length(sink_dim))
        }
        name <- sprintf("%09d", length(sink@dump))
        stopifnot(nchar(name) == 9L)
        assign(name, Rle(x), envir=sink@dump)
    }
)

setAs("RleRealizationSink", "RleArraySeed",
    function(from)
    {
        if (length(from@dump) == 0L) {
            rle <- Rle(get(from@type)(0))
        } else {
            rle <- do.call("c", unname(as.list(from@dump, sorted=TRUE)))
        }
        RleArraySeed(rle, from@dim, from@dimnames)
    }
)

setAs("RleRealizationSink", "RleArray",
    function(from) RleArray(as(from, "RleArraySeed"))
)

setAs("RleRealizationSink", "DelayedArray",
    function(from) RleArray(as(from, "RleArraySeed"))
)

.as_RleArray <- function(from)
{
    sink <- RleRealizationSink(dim(from), dimnames(from), type(from))
    write_to_sink(from, sink)
    as(sink, "RleArray")
}

setAs("ANY", "RleArray", .as_RleArray)

### Automatic coercion method from DelayedArray to RleArray silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "RleArray", .as_RleArray)

