### =========================================================================
### Read/write blocks of array data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_block() and write_block()
###

### Must return an ordinay array and propagate the dimnames.
setGeneric("read_block", signature="x",
    function(x, viewport)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)))
        ans <- standardGeneric("read_block")
        check_returned_array(ans, dim(viewport), "read_block", class(x))
    }
)

### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_block", signature="x",
    function(x, viewport, block)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is.array(block),
                  identical(dim(block), dim(viewport)))
        standardGeneric("write_block")
    }
)

### Work on any object 'x' that supports extract_array() e.g. an ordinary
### array or a DelayedArray, HDF5ArraySeed, or DelayedOp object, etc...
### Propagate the dimnames.
setMethod("read_block", "ANY",
    function(x, viewport)
    {
        ## We use extract_array_by_Nindex() instead of extract_array() to
        ## propagate the dimnames.
        Nindex <- makeNindexFromArrayViewport(viewport)
        extract_array_by_Nindex(x, Nindex)
    }
)

setMethod("write_block", "ANY",
    function(x, viewport, block)
    {
        Nindex <- makeNindexFromArrayViewport(viewport)
        replace_by_Nindex(x, Nindex, block)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block() and write_sparse_block()
###

### Like extract_sparse_array(), which it is based on, read_sparse_block()
### should be called only on an array-like object 'x' for which 'is_sparse(x)'
### is TRUE. For the sake of efficiency, this is NOT checked and is the
### responsibility of the user. See SparseArraySeed-class.R
### Must return a SparseArraySeed object.
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
        ans
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
### more efficient (can be 2x to 3x faster on big SparseArraySeed objects
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

### 'sparse_block' must be a SparseArraySeed object.
### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_sparse_block", signature="x",
    function(x, viewport, sparse_block)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is(sparse_block, "SparseArraySeed"),
                  identical(dim(sparse_block), dim(viewport)))
        standardGeneric("write_sparse_block")
    }
)

setMethod("write_sparse_block", "ANY",
    function(x, viewport, sparse_block)
    {
        block <- sparse2dense(sparse_block)
        write_block(x, viewport, block)
    }
)

