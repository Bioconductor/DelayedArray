### =========================================================================
### RealizationSink objects
### -------------------------------------------------------------------------


### Virtual class with no slots. Intended to be extended to support specific
### realization backends. Concrete subclasses must implement the "sink
### contract", that is:
###   1) A constructor function where the first 3 arguments are 'dim',
###      'dimnames', and 'type', in that order. Optionally it can have
###      the 'as.sparse' argument, in which case this **must** be the 4th
###      argument. It can have any additional argument.
###   2) A dim(), dimnames(), and type() method.
###   3) A write_block() method. It must return the modified array-like
###      object 'sink'.
###   4) A close() method (optional).
###   5) Coercion to DelayedArray.
###
### Examples of RealizationSink concrete subclasses: arrayRealizationSink
### (see below), RleRealizationSink (see RleArray-class.R),
### HDF5RealizationSink and TENxRealizationSink (see HDF5Array package).

setClass("RealizationSink", representation("VIRTUAL"))

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

setMethod("write_block", "arrayRealizationSink",
    function(sink, viewport, block)
    {
        result <- .get_arrayRealizationSink_result(sink)
        result <- write_block(result, viewport, block)
        .set_arrayRealizationSink_result(sink, result)
    }
)

setAs("arrayRealizationSink", "DelayedArray",
    function(from) DelayedArray(.get_arrayRealizationSink_result(from))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get/set the "realization backend" for the current session
###

.auto_realization_backend_envir <- new.env(parent=emptyenv())

getAutoRealizationBackend <- function()
{
    S4Arrays:::get_user_option("auto.realization.backend")
}

getRealizationBackend <- function()
{
    .Defunct("getAutoRealizationBackend")
    getAutoRealizationBackend()
}

.REGISTERED_REALIZATION_BACKENDS <- data.frame(
    BACKEND=c("RleArray", "HDF5Array", "TENxMatrix"),
    package=c("DelayedArray", "HDF5Array", "HDF5Array"),
    realization_sink_class=c("RleRealizationSink",
                             "HDF5RealizationSink",
                             "TENxRealizationSink"),
    stringsAsFactors=FALSE
)

registeredRealizationBackends <- function()
{
    ans <- .REGISTERED_REALIZATION_BACKENDS[ , c("BACKEND", "package")]
    backend <- getAutoRealizationBackend()
    Lcol <- ifelse(ans[ , "BACKEND"] %in% backend, "->", "")
    Rcol <- ifelse(ans[ , "BACKEND"] %in% backend, "<-", "")
    cbind(data.frame(` `=Lcol, check.names=FALSE),
          ans,
          data.frame(` `=Rcol, check.names=FALSE))
}

supportedRealizationBackends <- function()
{
    .Deprecated("registeredRealizationBackends")
    registeredRealizationBackends()
}

### NOT exported.
load_BACKEND_package <- function(BACKEND)
{
    if (!isSingleString(BACKEND))
        stop(wmsg("'BACKEND' must be a single string or NULL"))
    backends <- .REGISTERED_REALIZATION_BACKENDS
    m <- match(BACKEND, backends[ , "BACKEND"])
    if (is.na(m))
        stop(wmsg("Realization backend ", BACKEND, " is not registered ",
                  "in the DelayedArray package. Please use ",
                  "registeredRealizationBackends() to get the list of ",
                  "realization backends that are currently registered ",
                  "in the package."))
    package <- backends[ , "package"][[m]]
    class_package <- attr(BACKEND, "package")
    if (is.null(class_package)) {
        attr(BACKEND, "package") <- package
    } else if (!identical(package, class_package)) {
        stop(wmsg("\"package\" attribute on supplied 'BACKEND' is ",
                  "inconsistent with package normally associated with ",
                  "this realization backend"))
    }
    S4Arrays:::load_package_gracefully(package, "using the ",
                                       BACKEND, " realization backend")
    stopifnot(getClass(BACKEND)@package == package)
}

.check_realization_sink_constructor <- function(realization_sink_constructor)
{
    stopifnot(is.function(realization_sink_constructor))
    ok <- identical(head(formalArgs(realization_sink_constructor), n=3L),
                    c("dim", "dimnames", "type"))
    if (!ok)
        stop(wmsg("the first 3 arguments of a RealizationSink constructor ",
                  "function must be 'dim', 'dimnames', and 'type', in ",
                  "that order"))
    ## Either 'realization_sink_constructor' has the 'as.sparse' argument,
    ## in which case it **must** be in 4th position, or it does not have it.
    m <- match("as.sparse", formalArgs(realization_sink_constructor))
    if (!(m %in% c(4L, NA_integer_)))
        stop(wmsg("RealizationSink constructor functions with an 'as.sparse' ",
                  "argument must have it in 4th position"))
}

.get_realization_sink_constructor <- function(BACKEND)
{
    backends <- .REGISTERED_REALIZATION_BACKENDS
    m <- match(BACKEND, backends[ , "BACKEND"])
    realization_sink_class <- backends[ , "realization_sink_class"][[m]]
    package <- backends[ , "package"][[m]]
    realization_sink_constructor <- get(realization_sink_class,
                                        envir=.getNamespace(package),
                                        inherits=FALSE)
    .check_realization_sink_constructor(realization_sink_constructor)
    realization_sink_constructor
}

setAutoRealizationBackend <- function(BACKEND=NULL)
{
    if (is.null(BACKEND)) {
        remove(list=ls(envir=.auto_realization_backend_envir),
               envir=.auto_realization_backend_envir)
    } else {
        load_BACKEND_package(BACKEND)
        auto_realization_sink_constructor <-
            .get_realization_sink_constructor(BACKEND)
        assign("AUTO_REALIZATION_SINK_CONSTRUCTOR",
               auto_realization_sink_constructor,
               envir=.auto_realization_backend_envir)
    }
    S4Arrays:::set_user_option("auto.realization.backend", BACKEND)
    return(invisible(NULL))
}

setRealizationBackend <- function(...)
{
    .Defunct("setAutoRealizationBackend")
    setAutoRealizationBackend(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Backend-agnostic RealizationSink constructor
###

.get_auto_realization_sink_constructor <- function()
{
    if (is.null(getAutoRealizationBackend()))
        return(arrayRealizationSink)
    auto_realization_sink_constructor <-
        try(get("AUTO_REALIZATION_SINK_CONSTRUCTOR",
                envir=.auto_realization_backend_envir),
                silent=TRUE)
    if (is(auto_realization_sink_constructor, "try-error"))
        stop(wmsg("This operation requires a \"realization backend\". ",
                  "Please see '?setAutoRealizationBackend' for how ",
                  "to set one."))
    auto_realization_sink_constructor
}

.AutoRealizationSink <- function(dim, dimnames=NULL, type="double",
                                 as.sparse=FALSE, ...)
{
    realization_sink_constructor <- .get_auto_realization_sink_constructor()
    args <- list(dim, dimnames, type)
    formal_args <- formalArgs(realization_sink_constructor)
    if (length(formal_args) >= 4L && formal_args[[4L]] == "as.sparse")
        args <- c(args, list(as.sparse=as.sparse))
    args <- c(args, list(...))
    do.call(realization_sink_constructor, args)
}

AutoRealizationSink <- function(dim, dimnames=NULL, type="double",
                                as.sparse=FALSE)
{
    .AutoRealizationSink(dim, dimnames, type, as.sparse)
}

### Not exported.
RealizationSink <- function(BACKEND, ...)
{
    OLD_BACKEND <- getAutoRealizationBackend()
    setAutoRealizationBackend(BACKEND)
    on.exit(setAutoRealizationBackend(OLD_BACKEND))
    .AutoRealizationSink(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sinkApply()
###
### Thin wrapper around gridReduce().
###
### Note that, despite its name, this is actually a convenience wrapper
### around gridReduce(), and **not** around gridApply(). However, we call it
### sinkApply() and make its interface look similar to the gridApply/blockApply
### interface because this seems more user-friendly.
### Finally note that an important difference with gridReduce() is that the
### first two arguments of callback function 'FUN' are expected to be 'sink'
### and 'viewport' (in that order) rather than 'viewport' and 'init'.

sinkApply <- function(sink, FUN, ..., grid=NULL, verbose=NA)
{
    if (!is(sink, "RealizationSink"))
        stop(wmsg("'sink' must be a RealizationSink derivative"))
    FUN <- match.fun(FUN)
    grid <- normarg_sink_grid(grid, sink)
    verbose <- normarg_verbose(verbose)

    ## Only purpose of this wrapper is to swap the order of the first
    ## two arguments.
    FUN_WRAPPER <- function(viewport, init, FUN, ...)
    {
        effective_grid <- effectiveGrid()
        current_block_id <- currentBlockId()
        set_grid_context(effective_grid, current_block_id, viewport)
        FUN(init, viewport, ...)
    }
    gridReduce(FUN_WRAPPER, grid, sink, FUN, ..., verbose=verbose)
}

