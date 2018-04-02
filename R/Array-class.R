### =========================================================================
### Array objects
### -------------------------------------------------------------------------


### A virtual class with no slots to be extended by concrete subclasses with
### an array-like semantic.
setClass("Array", representation("VIRTUAL"))

### Even though prod() always returns a double, it seems that the length()
### primitive function takes care of turning this double into an integer if
### it's <= .Machine$integer.max
setMethod("length", "Array", function(x) prod(dim(x)))

setMethod("isEmpty", "Array", function(x) any(dim(x) == 0L))

### 'subscripts' is assumed to be an integer vector parallel to 'dim(x)' and
### with no out-of-bounds subscripts (i.e. 'all(subscripts >= 1)' and
### 'all(subscripts <= dim(x))').
### NOT exported for now but should probably be at some point (like
### S4Vectors::getListElement() is).
setGeneric("getArrayElement", signature="x",
    function(x, subscripts) standardGeneric("getArrayElement")
)

### Return an integer vector parallel to 'dim' and guaranteed to contain no
### out-of-bounds subscripts.
.from_linear_to_multi_subscript <- function(i, dim)
{
    stopifnot(isSingleInteger(i))
    if (i < 1L || i > prod(dim))
        stop("subscript is out of bounds")
    i <- i - 1L
    subscripts <- integer(length(dim))
    for (along in seq_along(dim)) {
        d <- dim[[along]]
        subscripts[[along]] <- offset <- i %% d
        i <- (i - offset) %/% d
    }
    subscripts + 1L
}

### Support multi-dimensional and linear subsetting.
### FIXME: Linear subsetting should support a single *numeric* subscript.
### FIXME: Multi-dimensional subsetting should support things like
###        x[[5, 15, 2]] and x[["E", 15, "b"]].
setMethod("[[", "Array",
    function(x, i, j, ...)
    {
        if (missing(x))
            stop("'x' is missing")
        Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
        nsubscript <- length(Nindex)
        x_dim <- dim(x)
        x_ndim <- length(x_dim)
        if (!(nsubscript == x_ndim || nsubscript == 1L))
            stop("incorrect number of subscripts")
        ok <- vapply(Nindex, isSingleInteger, logical(1), USE.NAMES=FALSE)
        if (!all(ok))
            stop(wmsg("each subscript must be a single integer ",
                      "when subsetting an ", class(x), " object with [["))
        if (nsubscript == x_ndim) {
            subscripts <- unlist(Nindex, use.names=FALSE)
            if (!(all(subscripts >= 1L) && all(subscripts <= x_dim)))
                stop("some subscripts are out of bounds")
        } else {
            ## Translate linear subsetting into multi-dimensional subsetting.
            subscripts <- .from_linear_to_multi_subscript(Nindex[[1L]], x_dim)
        }
        getArrayElement(x, subscripts)
    }
)

