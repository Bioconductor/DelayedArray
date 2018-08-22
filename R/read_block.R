### =========================================================================
### Read/write blocks of array data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block() and write_sparse_block()
###

### Must return a SparseData object.
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
        stopifnot(is(ans, "SparseData"),
                  identical(dim(ans), dim(viewport)))
        ans
    }
)

### 'sparse_block' must be a SparseData object.
### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_sparse_block", signature="x",
    function(x, viewport, sparse_block)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is(sparse_block, "SparseData"),
                  identical(dim(sparse_block), dim(viewport)))
        standardGeneric("write_sparse_block")
    }
)


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

