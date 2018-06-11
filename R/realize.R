### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RealizationSink objects
###

### Virtual class with no slots. Intended to be extended by implementations
### of DelayedArray backends. Concrete subclasses must implement:
###   1) A constructor function that takes argument 'dim', 'dimnames', and
###      'type'.
###   2) "dim" and "dimnames" methods.
###   3) A "write_block_to_sink" method.
###   4) A "close" method (optional).
###   5) Coercion to DelayedArray.
### See the arrayRealizationSink class below, or the RleRealizationSink class
### in RleArray-class.R, or the HDF5RealizationSink class in the HDF5Array
### package for examples of concrete RealizationSink subclasses.
setClass("RealizationSink", representation("VIRTUAL"))

### 'block', 'sink', and 'viewport' are expected to be an ordinary array, a
### RealizationSink, and a Viewport object, respectively. They must satisfy:
###
###    stopifnot(identical(dim(sink), refdim(viewport)),
###              identical(dim(block), dim(viewport)))
###
### Just to be safe, methods should perform this sanity check.
setGeneric("write_block_to_sink", signature="sink",
    function(block, sink, viewport) standardGeneric("write_block_to_sink")
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

.set_arrayRealizationSink_result <- function(sink, result)
{
    assign("result", result, envir=sink@result_envir)
}

setMethod("dim", "arrayRealizationSink",
    function(x) dim(.get_arrayRealizationSink_result(x))
)

arrayRealizationSink <- function(dim, dimnames=NULL, type="double")
{
    result <- array(get(type)(0), dim=dim, dimnames=dimnames)
    result_envir <- new.env(parent=emptyenv())
    sink <- new("arrayRealizationSink", result_envir=result_envir)
    .set_arrayRealizationSink_result(sink, result)
    sink
}

setMethod("write_block_to_sink", "arrayRealizationSink",
    function(block, sink, viewport)
    {
        stopifnot(identical(dim(sink), refdim(viewport)),
                  identical(dim(block), dim(viewport)))
        result <- .get_arrayRealizationSink_result(sink)
        result <- replace_block(result, viewport, block)
        .set_arrayRealizationSink_result(sink, result)
    }
)

setAs("arrayRealizationSink", "DelayedArray",
    function(from) DelayedArray(.get_arrayRealizationSink_result(from))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get/set the "realization backend" for the current session
###

.realization_backend_envir <- new.env(parent=emptyenv())

getRealizationBackend <- function()
{
    BACKEND <- try(get("BACKEND", envir=.realization_backend_envir),
                   silent=TRUE)
    if (is(BACKEND, "try-error"))
        return(NULL)
    BACKEND
}

.SUPPORTED_REALIZATION_BACKENDS <- data.frame(
    BACKEND=c("RleArray", "HDF5Array"),
    package=c("DelayedArray", "HDF5Array"),
    realization_sink_class=c("RleRealizationSink", "HDF5RealizationSink"),
    stringsAsFactors=FALSE
)

supportedRealizationBackends <- function()
{
    ans <- .SUPPORTED_REALIZATION_BACKENDS[ , c("BACKEND", "package")]
    backend <- getRealizationBackend()
    Lcol <- ifelse(ans[ , "BACKEND"] %in% backend, "->", "")
    Rcol <- ifelse(ans[ , "BACKEND"] %in% backend, "<-", "")
    cbind(data.frame(` `=Lcol, check.names=FALSE),
          ans,
          data.frame(` `=Rcol, check.names=FALSE))
}

.load_BACKEND_package <- function(BACKEND)
{
    if (!isSingleString(BACKEND))
        stop(wmsg("'BACKEND' must be a single string or NULL"))
    backends <- .SUPPORTED_REALIZATION_BACKENDS
    m <- match(BACKEND, backends[ , "BACKEND"])
    if (is.na(m))
        stop(wmsg("\"", BACKEND, "\" is not a supported backend. Please ",
                  "use supportedRealizationBackends() to get the list of ",
                  "supported \"realization backends\"."))
    package <- backends[ , "package"][[m]]
    class_package <- attr(BACKEND, "package")
    if (is.null(class_package)) {
        attr(BACKEND, "package") <- package
    } else if (!identical(package, class_package)) {
        stop(wmsg("\"package\" attribute on supplied 'BACKEND' is ",
                  "inconsistent with package normally associated with ",
                  "this backend"))
    }
    library(package, character.only=TRUE)
    stopifnot(getClass(BACKEND)@package == package)
}

.get_REALIZATION_SINK_CONSTRUCTOR <- function(BACKEND)
{
    backends <- .SUPPORTED_REALIZATION_BACKENDS
    m <- match(BACKEND, backends[ , "BACKEND"])
    realization_sink_class <- backends[ , "realization_sink_class"][[m]]
    package <- backends[ , "package"][[m]]
    REALIZATION_SINK_CONSTRUCTOR <- get(realization_sink_class,
                                        envir=.getNamespace(package),
                                        inherits=FALSE)
    stopifnot(is.function(REALIZATION_SINK_CONSTRUCTOR))
    stopifnot(identical(head(formalArgs(REALIZATION_SINK_CONSTRUCTOR), n=3L),
                        c("dim", "dimnames", "type")))
    REALIZATION_SINK_CONSTRUCTOR
}

setRealizationBackend <- function(BACKEND=NULL)
{
    if (is.null(BACKEND)) {
        remove(list=ls(envir=.realization_backend_envir),
               envir=.realization_backend_envir)
        return(invisible(NULL))
    }
    .load_BACKEND_package(BACKEND)
    REALIZATION_SINK_CONSTRUCTOR <- .get_REALIZATION_SINK_CONSTRUCTOR(BACKEND)
    assign("BACKEND", BACKEND,
           envir=.realization_backend_envir)
    assign("REALIZATION_SINK_CONSTRUCTOR", REALIZATION_SINK_CONSTRUCTOR,
           envir=.realization_backend_envir)
    return(invisible(NULL))
}

.get_realization_sink_constructor <- function()
{
    if (is.null(getRealizationBackend()))
        return(arrayRealizationSink)
    REALIZATION_SINK_CONSTRUCTOR <- try(get("REALIZATION_SINK_CONSTRUCTOR",
                                            envir=.realization_backend_envir),
                                        silent=TRUE)
    if (is(REALIZATION_SINK_CONSTRUCTOR, "try-error"))
        stop(wmsg("This operation requires a \"realization backend\". ",
                  "Please see '?setRealizationBackend' for how to set one."))
    REALIZATION_SINK_CONSTRUCTOR
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

setGeneric("realize", function(x, ...) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x, BACKEND=getRealizationBackend())
    {
        x <- DelayedArray(x)
        if (is.null(BACKEND))
            return(DelayedArray(as.array(x)))
        .load_BACKEND_package(BACKEND)
        ans <- as(x, BACKEND)
        ## Temporarily needed because coercion to HDF5Array currently drops
        ## the dimnames. See R/writeHDF5Array.R in the HDF5Array package for
        ## more information about this.
        ## TODO: Remove line below when this is addressed.
        set_dimnames(ans, dimnames(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RealizationSink constructor
###

RealizationSink <- function(dim, dimnames=NULL, type="double")
{
    .get_realization_sink_constructor()(dim, dimnames, type)
}

