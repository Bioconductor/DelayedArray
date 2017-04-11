### =========================================================================
### RealizationSink objects
### -------------------------------------------------------------------------
###
### Virtual class with no slots. Intended to be extended by implementations
### of DelayedArray backends. Concrete subclasses must implement:
###   1) A constructor function that takes argument 'dim', 'dimnames', and
###      'type'.
###   2) A "write_to_sink" method that works on an ordinary array.
###   3) A "close" method (optional).
###   4) Coercion to DelayedArray.
### See the arrayRealizationSink class below or the HDF5RealizationSink class
### in the HDF5Array package for examples of concrete RealizationSink
### subclasses.
###

setClass("RealizationSink", representation("VIRTUAL"))

### 'x' and 'sink' must have the same number of dimensions.
### 'offsets' must be NULL or an integer vector with 1 offset per dimension
### in 'x' (or in 'sink').
### A default "write_to_sink" method is defined in DelayedArray-class.R.
setGeneric("write_to_sink", signature=c("x", "sink"),
    function(x, sink, offsets=NULL) standardGeneric("write_to_sink")
)

setMethod("write_to_sink", c("DelayedArray", "RealizationSink"),
    function(x, sink, offsets=NULL)
    {
        if (!is.null(offsets))
            stop(wmsg("'offsets' must be NULL when the object to write ",
                      "is a DelayedArray object"))
        ## Semantically equivalent to 'write_to_sink(as.array(x), sink)'
        ## but uses block-processing so the full DelayedArray object is
        ## not realized at once in memory. Instead the object is first
        ## split into blocks and the blocks are realized and written to
        ## disk one at a time.
        block_APPLY(x, identity, sink=sink)
    }
)

setMethod("write_to_sink", c("ANY", "RealizationSink"),
    function(x, sink, offsets=NULL)
    {
        x <- as(x, "DelayedArray")
        write_to_sink(x, sink, offsets=offsets)
    }
)

setGeneric("close")

### The default "close" method for RealizationSink objects is a no-op.
setMethod("close", "RealizationSink", function(con) invisible(NULL))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### arrayRealizationSink objects
###
### The arrayRealizationSink class is a concrete RealizationSink subclass that
### implements an in-memory realization sink.
###

setClass("arrayRealizationSink",
    contains="RealizationSink",
    representation(
        result_envir="environment"
    )
)

.get_arrayRealizationSink_result <- function(sink)
{
    get("result", envir=sink@result_envir)
}

arrayRealizationSink <- function(dim, dimnames=NULL, type="double")
{
    result <- array(get(type)(0), dim=dim, dimnames=dimnames)
    result_envir <- new.env(parent=emptyenv())
    assign("result", result, envir=result_envir)
    new("arrayRealizationSink", result_envir=result_envir)
}

setMethod("write_to_sink", c("array", "arrayRealizationSink"),
    function(x, sink, offsets=NULL)
    {
        x_dim <- dim(x)
        result <- .get_arrayRealizationSink_result(sink)
        sink_dim <- dim(result)
        stopifnot(length(x_dim) == length(sink_dim))
        if (is.null(offsets)) {
            stopifnot(identical(x_dim, sink_dim))
            result[] <- x
        } else {
            block_ranges <- IRanges(offsets, width=x_dim)
            subscripts <- make_subscripts_from_block_ranges(
                              block_ranges, sink_dim,
                              expand.RangeNSBS=TRUE)
            result <- do.call(`[<-`, c(list(result), subscripts, list(value=x)))
        }
        assign("result", result, envir=sink@result_envir)
    }
)

setAs("arrayRealizationSink", "DelayedArray",
    function(from) DelayedArray(.get_arrayRealizationSink_result(from))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RealizationSink constructor
###

RealizationSink <- function(dim, dimnames=NULL, type="double")
{
    if (is.null(getRealizationBackend())) {
        REALIZATION_SINK_CONSTRUCTOR <- arrayRealizationSink
    } else {
        REALIZATION_SINK_CONSTRUCTOR <- get_realization_sink_constructor()
    }
    REALIZATION_SINK_CONSTRUCTOR(dim, dimnames, type)
}

