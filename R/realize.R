### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set/get the "realize() backend"
###

.SUPPORTED_REALIZE_BACKENDS <- data.frame(
    class="HDF5Array",
    package="HDF5Array",
    dump_constructor="HDF5ArrayDump",
    stringsAsFactors=FALSE
)

supportedRealizeBackends <- function()
    .SUPPORTED_REALIZE_BACKENDS[ , c("class", "package")]

.realize_backend_envir <- new.env(parent=emptyenv())

setRealizeBackend <- function(class=NULL)
{
    if (is.null(class)) {
        remove(list=ls(envir=.realize_backend_envir),
               envir=.realize_backend_envir)
        return(invisible(NULL))
    }
    if (!isSingleString(class))
        stop(wmsg("'class' must be a single string or NULL"))
    backends <- .SUPPORTED_REALIZE_BACKENDS
    m <- match(class, backends[ , "class"])
    if (is.na(m))
        stop(wmsg("\"", class, "\" is not a supported backend. Please use ",
                  "supportedRealizeBackends() to see the list of supported ",
                  "\"realize() backends\"."))
    package <- backends[ , "package"][[m]]
    dump_constructor <- backends[ , "dump_constructor"][[m]]
    class_package <- attr(class, "package")
    if (is.null(class_package)) {
        attr(class, "package") <- package
    } else if (!identical(package, class_package)) {
        stop(wmsg("\"package\" attribute on supplied 'class' is inconsistent ",
                  "with package normally associated with this backend"))
    }
    library(package, character.only=TRUE)
    stopifnot(getClass(class)@package == package)
    DUMP_CONSTRUCTOR <- get(dump_constructor, envir=.getNamespace(package),
                            inherits=FALSE)
    stopifnot(is.function(DUMP_CONSTRUCTOR))
    stopifnot(identical(formalArgs(DUMP_CONSTRUCTOR),
                        c("dim", "dimnames", "type")))
    assign("class", class, envir=.realize_backend_envir)
    assign("DUMP_CONSTRUCTOR", DUMP_CONSTRUCTOR, envir=.realize_backend_envir)
    return(invisible(NULL))
}

getRealizeBackend <- function()
{
    class <- try(get("class", envir=.realize_backend_envir), silent=TRUE)
    if (is(class, "try-error"))
        stop(wmsg("No \"realize() backend\" set yet. ",
                  "Please see '?setRealizeBackend' for how to set one."))
    class
}

### For internal use only (NOT exported).
get_DUMP_CONSTRUCTOR <- function()
{
    DUMP_CONSTRUCTOR <- try(get("DUMP_CONSTRUCTOR",
                                envir=.realize_backend_envir), silent=TRUE)
    if (is(DUMP_CONSTRUCTOR, "try-error"))
        stop(wmsg("This operation requires a \"realize() backend\". ",
                  "Please see '?setRealizeBackend' for how to set one."))
    DUMP_CONSTRUCTOR
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

realize <- function(x) as(DelayedArray(x), getRealizeBackend())

