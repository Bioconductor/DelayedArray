### =========================================================================
### SparseArraySeed objects
### -------------------------------------------------------------------------


setClass("SparseArraySeed",
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
### - Getters: dim(), length(), aind(), nzdata(), sparsity()
### - dense2sparse(), sparse2dense()
### - Based on sparse2dense(): extract_array(), as.array(), as.matrix()
### - Based on dense2sparse(): coercion to SparseArraySeed
### - Back and forth coercion between SparseArraySeed and dgCMatrix


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

.validate_SparseArraySeed <- function(x)
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

setValidity2("SparseArraySeed", .validate_SparseArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setGeneric("aind", function(x) standardGeneric("aind"))
setMethod("aind", "SparseArraySeed", function(x) x@aind)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "SparseArraySeed", function(x) x@nzdata)

setGeneric("sparsity", function(x) standardGeneric("sparsity"))
setMethod("sparsity", "SparseArraySeed",
    function(x) { 1 - length(nzdata(x)) / length(x) }
)


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

SparseArraySeed <- function(dim, aind=NULL, nzdata=NULL, check=TRUE)
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
    new2("SparseArraySeed", dim=dim, aind=aind, nzdata=nzdata, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object that supports 1D-style subsetting
### by a matrix like one returned by base::arrayInd(), that is, by a
### matrix where each row is an n-uplet representing an array index.
### Note that DelayedArray objects don't support this kind of subsetting
### yet so dense2sparse() doesn't work on them.
### Return a SparseArraySeed object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    aind <- which(x != 0L, arr.ind=TRUE)
    SparseArraySeed(x_dim, aind, x[aind], check=FALSE)
}

### 'sas' must be a SparseArraySeed object.
### Return an ordinary array.
sparse2dense <- function(sas)
{
    if (!is(sas, "SparseArraySeed"))
        stop(wmsg("'sas' must be a SparseArraySeed object"))
    ans <- array(0L, dim=dim(sas))
    ans[aind(sas)] <- nzdata(sas)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isSparse() and extract_sparse_array()
###

### isSparse() detects **structural** sparsity which is a qualitative
### property of array-like object 'x'. So it doesn't look at the data in 'x'.
### It is NOT about quantitative sparsity measured by sparsity().
setGeneric("isSparse", function(x) standardGeneric("isSparse"))

### By default, nothing is considered sparse.
setMethod("isSparse", "ANY", function(x) FALSE)

setMethod("isSparse", "SparseArraySeed", function(x) TRUE)

### Similar to extract_array() except that:
###   (1) The extracted array data must be returned in a SparseArraySeed
###       object.
###   (2) It should be called only on an array-like object 'x' for which
###       'isSparse(x)' is TRUE.
###   (3) The subscripts in 'index' should NOT contain duplicates.
### IMPORTANT NOTE: For the sake of efficiency, (2) and (3) are NOT checked
### and are the responsibility of the user. We'll refer to (2) and (3) as
### the "extract_sparse_array() Terms of Use".
setGeneric("extract_sparse_array",
    function(x, index)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("first argument to extract_sparse_array() ",
                      "must have dimensions"))
        ans <- standardGeneric("extract_sparse_array")
        expected_dim <- get_Nindex_lengths(index, x_dim)
        ## TODO: Display a more user/developper-friendly error by
        ## doing something like the extract_array() generic where
        ## check_returned_array() is used to display a long and
        ## detailed error message.
        stopifnot(is(ans, "SparseArraySeed"),
                  identical(dim(ans), expected_dim))
        ans
    }
)

### IMPORTANT NOTE: The returned SparseArraySeed object is guaranteed to be
### **correct** ONLY if the subscripts in 'index' do NOT contain duplicates!
### If they contain duplicates, the correct SparseArraySeed object to return
### should contain repeated nonzero data. However, in order to keep it as
### efficient as possible, the code below does NOT repeat the nonzero data
### that corresponds to duplicates subscripts. It does not check for
### duplicates in 'index' either because this check could have a
### non-negligible cost.
### All this is OK because .extract_sparse_array_from_SparseArraySeed()
### should always be used in a context where 'index' does NOT contain
### duplicates. The only situation where 'index' CAN contain duplicates
### is when .extract_sparse_array_from_SparseArraySeed() is called by
### .extract_array_from_SparseArraySeed(), in which case the missing
### nonzero data is added later.
.extract_sparse_array_from_SparseArraySeed <- function(x, index)
{
    stopifnot(is(x, "SparseArraySeed"))
    ans_dim <- get_Nindex_lengths(index, dim(x))
    x_aind <- x@aind
    for (along in seq_along(ans_dim)) {
        i <- index[[along]]
        if (is.null(i))
            next
        x_aind[ , along] <- match(x_aind[ , along], i)
    }
    keep_idx <- which(!rowAnyNAs(x_aind))
    ans_aind <- x_aind[keep_idx, , drop=FALSE]
    ans_nzdata <- x@nzdata[keep_idx]
    SparseArraySeed(ans_dim, ans_aind, ans_nzdata, check=FALSE)
}

setMethod("extract_sparse_array", "SparseArraySeed",
    .extract_sparse_array_from_SparseArraySeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_SparseArraySeed <- function(x, index)
{
    sas0 <- .extract_sparse_array_from_SparseArraySeed(x, index)
    ## If the subscripts in 'index' contain duplicates, 'sas0' is
    ## "incomplete" in the sense that it does not contain the nonzero data
    ## that should have been repeated according to the duplicates in the
    ## subscripts (see IMPORTANT NOTE above).
    ans0 <- sparse2dense(sas0)
    ## We "complete" 'ans0' by repeating the nonzero data according to the
    ## duplicates present in the subscripts. Note that this is easy and cheap
    ## to do now because 'ans0' uses a dense representation (it's an ordinary
    ## array). This would be much harder to do **natively** on the
    ## SparseArraySeed form (i.e. without converting first to dense then back
    ## to sparse in the process).
    sm_index <- lapply(index,
        function(i) {
            if (is.null(i))
                return(NULL)
            sm <- match(i, i)
            if (isSequence(sm))
                return(NULL)
            sm
        })
    if (all(S4Vectors:::sapply_isNULL(sm_index)))
        return(ans0)
    subset_by_Nindex(ans0, sm_index)
}

setMethod("extract_array", "SparseArraySeed",
    .extract_array_from_SparseArraySeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_sparse_block_from_SparseArraySeed()
###

### NOT exported but used in the HDF5Array package.
### Return a SparseArraySeed objects.
### TODO: Make this the "read_sparse_block" method for SparseArraySeed
### objects.
read_sparse_block_from_SparseArraySeed <- function(x, viewport)
{
    stopifnot(is(x, "SparseArraySeed"))
    taind <- t(x@aind)
    keep_idx <- which(colAlls(taind >= start(viewport) &
                              taind <= end(viewport)))
    taind <- taind[ , keep_idx, drop=FALSE]
    offsets <- start(viewport) - 1L
    x0_aind <- t(taind - offsets)
    x0_nzdata <- x@nzdata[keep_idx]
    BiocGenerics:::replaceSlots(x, dim=dim(viewport),
                                   aind=x0_aind,
                                   nzdata=x0_nzdata,
                                   check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to/from SparseArraySeed
###

### S3/S4 combo for as.array.SparseArraySeed
as.array.SparseArraySeed <- function(x, ...) sparse2dense(x)
setMethod("as.array", "SparseArraySeed", as.array.SparseArraySeed)

.from_SparseArraySeed_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    sparse2dense(x)
}

### S3/S4 combo for as.matrix.SparseArraySeed
as.matrix.SparseArraySeed <-
    function(x, ...) .from_SparseArraySeed_to_matrix(x, ...)
setMethod("as.matrix", "SparseArraySeed", .from_SparseArraySeed_to_matrix)

### Doesn't work on DelayedArray objects at the moment. See dense2sparse()
### above.
setAs("ANY", "SparseArraySeed", function(from) dense2sparse(from))

### Going back and forth between SparseArraySeed and dgCMatrix:

.from_dgCMatrix_to_SparseArraySeed <- function(from)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    aind <- cbind(i, j, deparse.level=0L)
    SparseArraySeed(dim(from), aind, from@x, check=FALSE)
}
setAs("dgCMatrix", "SparseArraySeed", .from_dgCMatrix_to_SparseArraySeed)

.from_SparseArraySeed_to_dgCMatrix <- function(from)
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
setAs("SparseArraySeed", "dgCMatrix", .from_SparseArraySeed_to_dgCMatrix)
setAs("SparseArraySeed", "sparseMatrix", .from_SparseArraySeed_to_dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###

.aperm.SparseArraySeed <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- normarg_perm(perm, dim(a))
    msg <- validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    BiocGenerics:::replaceSlots(a, dim=a_dim[perm],
                                   aind=a@aind[ , perm, drop=FALSE],
                                   check=FALSE)
}

### S3/S4 combo for aperm.SparseArraySeed
aperm.SparseArraySeed <-
    function(a, perm, ...) .aperm.SparseArraySeed(a, perm, ...)
setMethod("aperm", "SparseArraySeed", aperm.SparseArraySeed)

