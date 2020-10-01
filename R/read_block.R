### =========================================================================
### Read array blocks
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_block()
###

### Reads a block of data from array-like object 'x'. The block is returned
### either as an ordinay array or a SparseArraySeed object.
### Must propagate the dimnames. Note that this is taken care of by the
### generic function below so individual methods should not try to propagate
### the dimnames.
### 'as.sparse' can be TRUE, FALSE, or NA. If FALSE, the block is returned
### as an ordinary (dense) array. If TRUE, it's returned as a SparseArraySeed
### object. Using NA is equivalent to using 'as.sparse=is_sparse(x)' and is
### the most efficient way to read a block. Maybe it should be made the
### default (right now the default is FALSE for backward compatibility with
### existing code).
### IMPORTANT: Methods should ignore the 'as.sparse' argument! (i.e. have
### it in their argument list but not do anything about it).
setGeneric("read_block", signature="x",
    function(x, viewport, as.sparse=FALSE)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is.logical(as.sparse),
                  length(as.sparse) == 1L)
        if (is_sparse(x)) {
            ## read_sparse_block() is already taking care of propagating
            ## the dimnames.
            ans <- read_sparse_block(x, viewport)
            if (isFALSE(as.sparse))
                ans <- sparse2dense(ans)
        } else {
            ans <- standardGeneric("read_block")
            check_returned_array(ans, dim(viewport), "read_block", class(x))
            ## Individual read_block() methods should not try to propagate
            ## the dimnames. We do that for them now.
            Nindex <- makeNindexFromArrayViewport(viewport)
            ans_dimnames <- subset_dimnames_by_Nindex(dimnames(x), Nindex)
            ans <- set_dimnames(ans, ans_dimnames)
            if (isTRUE(as.sparse))
                ans <- dense2sparse(ans)
        }
        ans
    }
)

### Work on any object 'x' that supports extract_array() e.g. an ordinary
### array, a sparseMatrix derivative from the Matrix package, a DelayedArray
### object, an HDF5ArraySeed object, a DelayedOp object, etc...
setMethod("read_block", "ANY",
    function(x, viewport)
    {
        Nindex <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        extract_array(x, Nindex)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block()
###

### Like extract_sparse_array(), which it is based on, read_sparse_block()
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
        ans_dimnames <- subset_dimnames_by_Nindex(dimnames(x), Nindex)
        set_dimnames(ans, ans_dimnames)
    }
)

setMethod("read_sparse_block", "ANY",
    function(x, viewport)
    {
        index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        extract_sparse_array(x, index)
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

