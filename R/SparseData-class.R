### =========================================================================
### SparseData objects
### -------------------------------------------------------------------------


### Do NOT extend Array! The length() of an Array derivative 'x' is
### defined as 'prod(dim(x))' whereas here we want it to be something else.
### Note that there are other objects with dimensions that have a length()
### that is not 'prod(dim(x))' e.g. data-frame-like objects (for which
### 'length(x)' is 'ncol(x)') and SummarizedExperiment derivatives (for
### which 'length(x)' is 'nrow(x)') so we are not creating a precedent.
### Terminology: Should we still consider these objects to be "array-like"
### or "matrix-like"? Or should these terms be used only for objects that
### have dimensions **and** have a length() defined as 'prod(dim(x))'?
setClass("SparseData",
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
### - Based on dense2sparse(): coercion to SparseData
### - Based on sparse2dense(): as.array()
### - Other coercions: back and forth between SparseData and dgCMatrix.


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

.validate_SparseData <- function(x)
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

setValidity2("SparseData", .validate_SparseData)


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

SparseData <- function(dim, aind=NULL, nzdata=NULL, check=TRUE)
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
        nzdata <- .normarg_nzdata(nzdata, nrow(aind))
    }
    new2("SparseData", dim=dim, aind=aind, nzdata=nzdata, check=check)
}


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
    SparseData(x_dim, aind, x[aind], check=FALSE)
}

### 'sparse_data' must be a SparseData object.
### Return an ordinary array.
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
    ans_aind <- cbind(x@i + 1L,
                      rep.int(seq_len(ncol(x)), diff(x@p)),
                      deparse.level=0L)
    SparseData(dim(x), ans_aind, x@x, check=FALSE)
}

