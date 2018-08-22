### =========================================================================
### SparseData objects
### -------------------------------------------------------------------------


### Do NOT extend Array! The length() of an Array derivative 'x' is
### defined as 'prod(dim(x))' whereas here we want it to be something else.
### Other objects with a length() that disagrees with the length of an
### array derivative: data-frame-like objects ('length(x)' is 'ncol(x)')
### and SummarizedExperiment derivatives ('length(x)' is 'nrow(x)').
setClass("SparseData",
    representation(
        dim="integer",   # This gives us dim() for free!
        aind="matrix",   # A numeric matrix like one returned by
                         # base::arrayInd(), that is, a matrix with
                         # 'length(dim)' columns where each row is an
                         # n-uplet representing an array index.
        nzdata="vector"  # A vector of length 'nrow(aind)' containing the
                         # nonzero data.
    )
)

### API:
### - Getters: dim(), length(), aind(), nonzeroes()
### - dense2sparse(), sparse2dense()
### - Based on dense2sparse(): coercion to SparseData
### - Based on sparse2dense(): as.array()
### - Other coercions: back and forth between SparseData and dgCMatrix.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

### In disagreement with the length() of Array derivatives!
setMethod("length", "SparseData", function(x) length(x@nzdata))

setGeneric("aind", function(x) standardGeneric("aind"))
setMethod("aind", "SparseData", function(x) x@aind)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "SparseData", function(x) x@nzdata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object.
### Return a SparseData object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    aind <- which(x != 0L, arr.ind=TRUE)
    new2("SparseData", dim=x_dim, aind=aind, nzdata=x[aind], check=FALSE)
}

### 'sparse_data' must be a SparseData object.
### Return an ordinary matrix.
sparse2dense <- function(sparse_data)
{
    if (!is(sparse_data, "SparseData"))
        stop(wmsg("'sparse_data' must be a SparseData object"))
    ans <- array(0L, dim=dim(sparse_data))
    ans[aind(sparse_data)] <- nzdata(sparse_data)
    ans
}

### NOT exported and currently unused.
dense2sparse.dgCMatrix <- function(x)
{
    stopifnot(is(x, "dgCMatrix"))
    data.frame(i=x@i + 1L, j=rep.int(seq_len(ncol(x)), diff(x@p)), data=x@x)
}

