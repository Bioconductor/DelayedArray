### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get/set the "realize() backend"
###

.SUPPORTED_REALIZE_BACKENDS <- data.frame(
    BACKEND="HDF5Array",
    package="HDF5Array",
    realization_sink_class="HDF5RealizationSink",
    stringsAsFactors=FALSE
)

supportedRealizeBackends <- function()
    .SUPPORTED_REALIZE_BACKENDS[ , c("BACKEND", "package")]

.realize_backend_envir <- new.env(parent=emptyenv())

getRealizeBackend <- function()
{
    BACKEND <- try(get("BACKEND", envir=.realize_backend_envir), silent=TRUE)
    if (is(BACKEND, "try-error"))
        return(NULL)
    BACKEND
}

.load_BACKEND_package <- function(BACKEND)
{
    if (!isSingleString(BACKEND))
        stop(wmsg("'BACKEND' must be a single string or NULL"))
    backends <- .SUPPORTED_REALIZE_BACKENDS
    m <- match(BACKEND, backends[ , "BACKEND"])
    if (is.na(m))
        stop(wmsg("\"", BACKEND, "\" is not a supported backend. Please use ",
                  "supportedRealizeBackends() to see the list of supported ",
                  "\"realize() backends\"."))
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
    backends <- .SUPPORTED_REALIZE_BACKENDS
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

setRealizeBackend <- function(BACKEND=NULL)
{
    if (is.null(BACKEND)) {
        remove(list=ls(envir=.realize_backend_envir),
               envir=.realize_backend_envir)
        return(invisible(NULL))
    }
    .load_BACKEND_package(BACKEND)
    REALIZATION_SINK_CONSTRUCTOR <- .get_REALIZATION_SINK_CONSTRUCTOR(BACKEND)
    assign("BACKEND", BACKEND,
           envir=.realize_backend_envir)
    assign("REALIZATION_SINK_CONSTRUCTOR", REALIZATION_SINK_CONSTRUCTOR,
           envir=.realize_backend_envir)
    return(invisible(NULL))
}

get_realization_sink_constructor <- function()
{
    REALIZATION_SINK_CONSTRUCTOR <- try(get("REALIZATION_SINK_CONSTRUCTOR",
                                            envir=.realize_backend_envir),
                                        silent=TRUE)
    if (is(REALIZATION_SINK_CONSTRUCTOR, "try-error"))
        stop(wmsg("This operation requires a \"realize() backend\". ",
                  "Please see '?setRealizeBackend' for how to set one."))
    REALIZATION_SINK_CONSTRUCTOR
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

setGeneric("realize", function(x, ...) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x, BACKEND=getRealizeBackend())
    {
        x <- DelayedArray(x)
        if (is.null(BACKEND))
            return(as.array(x))
        .load_BACKEND_package(BACKEND)
        ans <- as(x, BACKEND)
        ## Temporarily needed because coercion to HDF5Array currently drops
        ## the dimnames. See R/writeHDF5Array.R in the HDF5Array package for
        ## more information about this.
        ## TODO: Remove line below when this is addressed.
        dimnames(ans) <- dimnames(x)
        ans
    }
)

