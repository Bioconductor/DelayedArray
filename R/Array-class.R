### =========================================================================
### Array objects
### -------------------------------------------------------------------------


### A virtual class with no slots to be extended by concrete subclasses with
### an array-like semantic.
setClass("Array", representation("VIRTUAL"))

### Note that some objects with dimensions (i.e. with a non-NULL dim()) can
### have a length() that is not 'prod(dim(x))' e.g. data-frame-like objects
### (for which 'length(x)' is 'ncol(x)') and SummarizedExperiment derivatives
### (for which 'length(x)' is 'nrow(x)').
### Terminology: Should we still consider that these objects are "array-like"
### or "matrix-like"? Or should these terms be used only for objects that
### have dimensions **and** have a length() defined as 'prod(dim(x))'?

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

### Support multi-dimensional and linear subsetting.
### TODO: Multi-dimensional subsetting should support things like
###       x[[5, 15, 2]] and x[["E", 15, "b"]].
### TODO: Linear subsetting should support a single *numeric* subscript.
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
            ## Multi-dimensional subsetting.
            subscripts <- unlist(Nindex, use.names=FALSE)
            if (!(all(subscripts >= 1L) && all(subscripts <= x_dim)))
                stop("some subscripts are out of bounds")
        } else {
            ## Linear subsetting.
            ## We turn this into a multi-dimensional subsetting by
            ## transforming the user-supplied linear index into an array
            ## (i.e. multi-dimensional) index.
            i <- Nindex[[1L]]
            if (i < 1L || i > prod(x_dim))
                stop("subscript is out of bounds")
            subscripts <- as.integer(arrayInd(i, x_dim))
        }
        getArrayElement(x, subscripts)
    }
)

### S3/S4 combo for t.Array
### t() will work out-of-the-box on any Array derivative that supports aperm().
t.Array <- function(x)
{
    if (length(dim(x)) != 2L)
        stop(wmsg("the ", class(x), " object to transpose ",
                  "must have exactly 2 dimensions"))
    aperm(x)
}
setMethod("t", "Array", t.Array)

### Any Array derivative will be coercible to a sparseMatrix subclass (e.g.
### dg[C|R]Matrix or lg[C|R]Matrix) as long as there is a method for coercing
### it to SparseArraySeed.
setAs("Array", "dgCMatrix",
    function(from) as(as(from, "SparseArraySeed"), "dgCMatrix")
)
setAs("Array", "dgRMatrix",
    function(from) as(as(from, "SparseArraySeed"), "dgRMatrix")
)
setAs("Array", "lgCMatrix",
    function(from) as(as(from, "SparseArraySeed"), "lgCMatrix")
)
setAs("Array", "lgRMatrix",
    function(from) as(as(from, "SparseArraySeed"), "lgRMatrix")
)
setAs("Array", "CsparseMatrix",
    function(from) as(as(from, "SparseArraySeed"), "CsparseMatrix")
)
setAs("Array", "RsparseMatrix",
    function(from) as(as(from, "SparseArraySeed"), "RsparseMatrix")
)
setAs("Array", "sparseMatrix",
    function(from) as(as(from, "SparseArraySeed"), "sparseMatrix")
)

