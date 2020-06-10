### =========================================================================
### Read/write blocks of array data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_block()
###

### Reads a block of data from array-like object 'x'. The block is returned
### either as an ordinay array or a SparseArraySeed object.
### Is expected to propagate the dimnames (mmh.. why is that a requirement?
### can't remember, maybe this could be relaxed).
### 'as.sparse' can be TRUE, FALSE, or NA. If FALSE, the block is returned
### as an ordinary (dense) array. If TRUE, it's returned as a SparseArraySeed
### object. Using NA is equivalent to using 'as.sparse=is_sparse(x)' and is
### the most efficient way to read a block. Maybe it should be made the
### default (right now the default is FALSE for backward compatibility with
### existing code).
### IMPORTANT: Methods should ignore the 'as.sparse' argument! (i.e. have
### it in their argument but not do anything about it).
setGeneric("read_block", signature="x",
    function(x, viewport, as.sparse=FALSE)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is.logical(as.sparse),
                  length(as.sparse) == 1L)
        if (is_sparse(x)) {
            ans <- read_sparse_block(x, viewport)
            if (isFALSE(as.sparse))
                ans <- sparse2dense(ans)
        } else {
            ans <- standardGeneric("read_block")
            check_returned_array(ans, dim(viewport), "read_block", class(x))
            if (isTRUE(as.sparse))
                ans <- dense2sparse(ans)
        }
        ans
    }
)

### Work on any object 'x' that supports extract_array() e.g. an ordinary
### array, a sparseMatrix derivative from the Matrix package, a DelayedArray
### object, an HDF5ArraySeed object, a DelayedOp object, etc...
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block()
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write_block()
###

### 'x' is typically expected to be a RealizationSink derivative but
### write_block() should also work on an ordinary array or other
### in-memory array- or matrix-like object like a sparseMatrix derivative
### from the Matrix package.
### Dispatch on first argument 'x' only for now for simplicity but we
### could change this to also dispatch on the third argument ('block')
### when the need arises.
### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_block", signature="x",
    function(x, viewport, block)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  identical(dim(block), dim(viewport)))
        standardGeneric("write_block")
    }
)

setMethod("write_block", "ANY",
    function(x, viewport, block)
    {
        if (is(block, "SparseArraySeed"))
            block <- sparse2dense(block)
        Nindex <- makeNindexFromArrayViewport(viewport)
        replace_by_Nindex(x, Nindex, block)
    }
)

