### =========================================================================
### Read/write blocks of array data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### Return a data.frame with 3 columns ("i", "j", "data") and 1 row per
### nonzero element in matrix-like object 'x'.
dense2sparse <- function(x)
{
    if (length(dim(x)) != 2L)
        stop(wmsg("'x' must be a matrix-like object"))
    aind <- which(x != 0L, arr.ind=TRUE)
    data.frame(i=aind[ , "row"], j=aind[ , "col"], data=x[aind],
               stringsAsFactors=FALSE)
}

### NOT exported and currently unused.
dense2sparse.dgCMatrix <- function(x)
{
    stopifnot(is(x, "dgCMatrix"))
    data.frame(i=x@i + 1L, j=rep.int(seq_len(ncol(x)), diff(x@p)), data=x@x)
}

### 'sparse_data' must be a data.frame as returned by dense2sparse().
### Return an ordinary matrix.
sparse2dense <- function(dim, sparse_data)
{
    if (!(is.numeric(dim) && length(dim) == 2L))
        stop(wmsg("'dim' must have length 2"))
    ans <- array(0L, dim=dim)
    ans[cbind(sparse_data$i, sparse_data$j)] <- sparse_data$data
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block() and write_sparse_block()
###
### Only supports matrix-like objects for now.
###

### Must return a data.frame as returned by dense2sparse().
setGeneric("read_sparse_block", signature="x",
    function(x, viewport)
    {
        x_dim <- dim(x)
        stopifnot(length(x_dim) == 2L,
                  is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), x_dim))
        standardGeneric("read_sparse_block")
    }
)

### 'sparse_block' must be a data.frame as returned by dense2sparse().
### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_sparse_block", signature="x",
    function(x, viewport, sparse_block)
    {
        x_dim <- dim(x)
        stopifnot(length(x_dim) == 2L,
                  is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), x_dim))
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

