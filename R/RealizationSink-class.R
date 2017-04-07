### =========================================================================
### Dump management utilities
### -------------------------------------------------------------------------
###

### Virtual class with no slots. Intended to be extended by implementations
### of DelayedArray on-disk backends. Concrete subclasses must implement:
###   1) A constructor function that takes argument 'dim', 'dimnames', and
###      'type'.
###   2) A "write_to_sink" method that works on an ordinary array.
###   3) A "close" method.
###   4) Coercion to DelayedArray.
### See HDF5RealizationSink class in the HDF5Array package for an example.
setClass("RealizationSink", representation("VIRTUAL"))

### A default "write_to_sink" method is defined in DelayedArray-class.R.
setGeneric("write_to_sink", signature=c("x", "sink"),
    function(x, sink, offsets=NULL) standardGeneric("write_to_sink")
)

setGeneric("close")

### The default "close" method for RealizationSink objects is a no-op.
setMethod("close", "RealizationSink", function(con) invisible(NULL))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RealizationSink constructor
###

.get_realization_sink_constructor <- function()
{
    REALIZATION_SINK_CONSTRUCTOR <- try(get("REALIZATION_SINK_CONSTRUCTOR",
                                            envir=.realize_backend_envir),
                                        silent=TRUE)
    if (is(REALIZATION_SINK_CONSTRUCTOR, "try-error"))
        stop(wmsg("This operation requires a \"realize() backend\". ",
                  "Please see '?setRealizeBackend' for how to set one."))
    REALIZATION_SINK_CONSTRUCTOR
}

RealizationSink <- function(dim, dimnames, type)
{
    .get_realization_sink_constructor()(dim, dimnames, type)
}

