### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


### Note that there are objects with dimensions that have a length()
### that is not 'prod(dim(x))' e.g. data-frame-like objects (for which
### 'length(x)' is 'ncol(x)') and SummarizedExperiment derivatives (for
### which 'length(x)' is 'nrow(x)').
### Terminology: Should we still consider these objects to be "array-like"
### or "matrix-like"? Or should these terms be used only for objects that
### have dimensions **and** have a length() defined as 'prod(dim(x))'?
setClass("SparseArray",
    contains="Array",
    representation(
        dim="integer",   # This gives us dim() for free!
        aind="matrix",   # An **integer** matrix like one returned by
                         # base::arrayInd(), that is, with 'length(dim)'
                         # columns and where each row is an n-uplet
                         # representing an "array index".
        nzdata="vector"  # A vector of length 'nrow(aind)' containing the
                         # nonzero data.
    )
)

### API:
### - Getters: dim(), length(), aind(), nonzeroes()
### - dense2sparse(), sparse2dense()
### - Based on dense2sparse(): coercion to SparseArray
### - Based on sparse2dense(): as.array(), as.matrix()
### - Other coercions: back and forth between SparseArray and dgCMatrix.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_aind_slot <- function(x)
{
    x_aind <- x@aind
    if (!(is.matrix(x_aind) && typeof(x_aind) == "integer"))
        return(wmsg2("'aind' slot must be an integer matrix"))
    x_dim <- x@dim
    if (ncol(x_aind) != length(x_dim))
        return(wmsg2("'aind' slot must be a matrix with ",
                     "one column per dimension"))
    for (along in seq_along(x_dim)) {
        notok <- S4Vectors:::anyMissingOrOutside(x_aind[ , along],
                                                 1L, x_dim[[along]])
        if (notok)
            return(wmsg2("'aind' slot must contain valid indices, ",
                         "that is, indices that are not NA and are ",
                         ">= 1 and <= their corresponding dimension"))
    }
    TRUE
}

.validate_nzdata_slot <- function(x)
{
    x_nzdata <- x@nzdata
    if (!(is.vector(x_nzdata) && length(x_nzdata) == nrow(x@aind)))
        return(wmsg2("'nzdata' slot must be a vector of length ",
                     "the number of rows in the 'aind' slot"))
    TRUE
}

.validate_SparseArray <- function(x)
{
    msg <- validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_aind_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzdata_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("SparseArray", .validate_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setGeneric("aind", function(x) standardGeneric("aind"))
setMethod("aind", "SparseArray", function(x) x@aind)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "SparseArray", function(x) x@nzdata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_nzdata <- function(nzdata, length.out)
{
    if (is.null(nzdata))
        stop(wmsg("'nzdata' cannot be NULL when 'aind' is not NULL"))
    if (!is.vector(nzdata))
        stop(wmsg("'nzdata' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzdata_len <- length(nzdata)
    if (nzdata_len == length.out)
        return(nzdata)
    if (nzdata_len > length.out && nzdata_len != 1L)
        stop(wmsg("'length(nzdata)' is greater than 'nrow(aind)'"))
    if (nzdata_len == 0L)
        stop(wmsg("'length(nzdata)' is 0 but 'nrow(aind)' is not"))
    if (length.out %% nzdata_len != 0L)
        warning(wmsg("'nrow(aind)' is not a multiple of 'length(nzdata)'"))
    rep(nzdata, length.out=length.out)
}

SparseArray <- function(dim, aind=NULL, nzdata=NULL, check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(aind)) {
        if (!is.null(nzdata))
            stop(wmsg("'nzdata' must be NULL when 'aind' is NULL"))
        aind <- matrix(integer(0), ncol=length(dim))
        nzdata <- integer(0)
    } else {
        if (!is.matrix(aind))
            stop(wmsg("'aind' must be a matrix"))
        if (storage.mode(aind) == "double")
            storage.mode(aind) <- "integer"
        if (!is.null(dimnames(aind)))
            dimnames(aind) <- NULL
        nzdata <- .normarg_nzdata(nzdata, nrow(aind))
    }
    new2("SparseArray", dim=dim, aind=aind, nzdata=nzdata, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object that supports 1D-style subsetting
### by a matrix like one returned by base::arrayInd(), that is, by a
### matrix where each row is an n-uplet representing an array index.
### Note that DelayedArray objects don't support this kind of subsetting
### yet so dense2sparse() doesn't work on them.
### Return a SparseArray object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    aind <- which(x != 0L, arr.ind=TRUE)
    SparseArray(x_dim, aind, x[aind], check=FALSE)
}

### 'sparse_array' must be a SparseArray object.
### Return an ordinary array.
sparse2dense <- function(sparse_array)
{
    if (!is(sparse_array, "SparseArray"))
        stop(wmsg("'sparse_array' must be a SparseArray object"))
    ans <- array(0L, dim=dim(sparse_array))
    ans[aind(sparse_array)] <- nzdata(sparse_array)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to/from SparseArray
###

### S3/S4 combo for as.array.SparseArray
as.array.SparseArray <- function(x, ...) sparse2dense(x)
setMethod("as.array", "SparseArray", as.array.SparseArray)

.from_SparseArray_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    sparse2dense(x)
}

### S3/S4 combo for as.matrix.SparseArray
as.matrix.SparseArray <- function(x, ...) .from_SparseArray_to_matrix(x, ...)
setMethod("as.matrix", "SparseArray", .from_SparseArray_to_matrix)

### Doesn't work on DelayedArray objects at the moment. See dense2sparse()
### above.
setAs("ANY", "SparseArray", function(from) dense2sparse(from))

### Going back and forth between SparseArray and dgCMatrix:

.from_dgCMatrix_to_SparseArray <- function(from)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    aind <- cbind(i, j, deparse.level=0L)
    SparseArray(dim(from), aind, from@x, check=FALSE)
}
setAs("dgCMatrix", "SparseArray", .from_dgCMatrix_to_SparseArray)

.from_SparseArray_to_dgCMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "must have exactly 2 dimensions"))
    i <- from@aind[ , 1L]
    j <- from@aind[ , 2L]
    x <- from@nzdata
    Matrix::sparseMatrix(i, j, x=x, dims=from_dim, dimnames=dimnames(from))
}
setAs("SparseArray", "dgCMatrix", .from_SparseArray_to_dgCMatrix)
setAs("SparseArray", "sparseMatrix", .from_SparseArray_to_dgCMatrix)

