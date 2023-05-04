### =========================================================================
### read_sparse_block()
### -------------------------------------------------------------------------
###

### Like OLD_extract_sparse_array(), which it is based on, read_sparse_block()
### should be called only on an array-like object 'x' for which 'is_sparse(x)'
### is TRUE. For the sake of efficiency, this is NOT checked and is the
### responsibility of the user. See SparseArraySeed-class.R
### Must return a SparseArraySeed object with the dimnames on 'x' propagated
### to it. Note that this is taken care of by the generic function below so
### individual methods should not try to propagate the dimnames.
setGeneric("read_sparse_block", signature="x",
    function(x, viewport)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)))
        ans <- standardGeneric("read_sparse_block")
        ## TODO: Display a more user/developper-friendly error by
        ## doing something like the read_block() generic where
        ## check_returned_array() is used to display a long and
        ## detailed error message.
        stopifnot(is(ans, "SparseArraySeed"),
                  identical(dim(ans), dim(viewport)))
        Nindex <- makeNindexFromArrayViewport(viewport)
        ans_dimnames <- S4Arrays:::subset_dimnames_by_Nindex(dimnames(x),
                                                             Nindex)
        S4Arrays:::set_dimnames(ans, ans_dimnames)
    }
)

setMethod("read_sparse_block", "ANY",
    function(x, viewport)
    {
        index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        OLD_extract_sparse_array(x, index)
    }
)

### The default "read_sparse_block" method above would work just fine on a
### SparseArraySeed object but we overwrite it with a method that is slightly
### more efficient (can be 2x to 3x faster on a big SparseArraySeed object
### with hundreds of thousands of nonzero elements).
.read_sparse_block_from_SparseArraySeed <- function(x, viewport)
{
    tnzindex <- t(x@nzindex)
    keep_idx <- which(colAlls(tnzindex >= start(viewport) &
                              tnzindex <= end(viewport)))
    tnzindex <- tnzindex[ , keep_idx, drop=FALSE]
    offsets <- start(viewport) - 1L
    x0_nzindex <- t(tnzindex - offsets)
    x0_nzdata <- x@nzdata[keep_idx]
    BiocGenerics:::replaceSlots(x, dim=dim(viewport),
                                   nzindex=x0_nzindex,
                                   nzdata=x0_nzdata,
                                   check=FALSE)
}

setMethod("read_sparse_block", "SparseArraySeed",
    .read_sparse_block_from_SparseArraySeed
)

