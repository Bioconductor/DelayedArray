### =========================================================================
### SparseArraySeed objects
### -------------------------------------------------------------------------


setClass("SparseArraySeed",
    contains="Array",
    representation(
        dim="integer",     # This gives us dim() for free!
        nzindex="matrix",  # M-index of the nonzero data.
        nzdata="vector",   # A vector (atomic or list) of length
                           # 'nrow(nzindex)' containing the nonzero data.
        dimnames="list"    # List with one list element per dimension. Each
                           # list element must be NULL or a character vector.
    ),
    prototype(
        dim=0L,
        nzindex=matrix(integer(0), ncol=1L),
        dimnames=list(NULL)
    )
)

### API:
### - Getters/setters: dim(), length(), nzindex(), nzdata(),
###                    dimnames(), dimnames<-()
### - sparsity().
### - dense2sparse(), sparse2dense().
### - Based on sparse2dense(): extract_array(), as.array(), as.matrix().
### - Based on dense2sparse(): coercion to SparseArraySeed.
### - Back and forth coercion between SparseArraySeed and dg[C|R]Matrix or
###   lg[C|R]Matrix objects from the Matrix package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_nzindex_slot <- function(x)
{
    x_nzindex <- x@nzindex
    if (!(is.matrix(x_nzindex) && typeof(x_nzindex) == "integer"))
        return("'nzindex' slot must be an integer matrix")
    x_dim <- x@dim
    if (ncol(x_nzindex) != length(x_dim))
        return(paste0("'nzindex' slot must be a matrix with ",
                      "one column per dimension"))
    for (along in seq_along(x_dim)) {
        not_ok <- S4Vectors:::anyMissingOrOutside(x_nzindex[ , along],
                                                  1L, x_dim[[along]])
        if (not_ok)
            return(paste0("'nzindex' slot must contain valid indices, ",
                          "that is, indices that are not NA and are ",
                          ">= 1 and <= their corresponding dimension"))
    }
    TRUE
}

.validate_nzdata_slot <- function(x)
{
    x_nzdata <- x@nzdata
    if (!(is.vector(x_nzdata) && length(x_nzdata) == nrow(x@nzindex)))
        return(paste0("'nzdata' slot must be a vector of length ",
                      "the number of rows in the 'nzindex' slot"))
    TRUE
}

.validate_SparseArraySeed <- function(x)
{
    msg <- validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzindex_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzdata_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dimnames_slot(x, x@dim)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("SparseArraySeed", .validate_SparseArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters
###

setGeneric("nzindex", function(x) standardGeneric("nzindex"))
setMethod("nzindex", "SparseArraySeed", function(x) x@nzindex)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "SparseArraySeed", function(x) x@nzdata)

setMethod("dimnames", "SparseArraySeed",
    function(x) simplify_NULL_dimnames(x@dimnames)
)

setReplaceMethod("dimnames", "SparseArraySeed",
    function(x, value)
    {
        x@dimnames <- normarg_dimnames(value, dim(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparsity()
###

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
        stop(wmsg("'nzdata' cannot be NULL when 'nzindex' is not NULL"))
    if (!is.vector(nzdata))
        stop(wmsg("'nzdata' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzdata_len <- length(nzdata)
    if (nzdata_len == length.out)
        return(nzdata)
    if (nzdata_len > length.out && nzdata_len != 1L)
        stop(wmsg("'length(nzdata)' is greater than 'nrow(nzindex)'"))
    if (nzdata_len == 0L)
        stop(wmsg("'length(nzdata)' is 0 but 'nrow(nzindex)' is not"))
    if (length.out %% nzdata_len != 0L)
        warning(wmsg("'nrow(nzindex)' is not a multiple of 'length(nzdata)'"))
    rep(nzdata, length.out=length.out)
}

SparseArraySeed <- function(dim, nzindex=NULL, nzdata=NULL, dimnames=NULL,
                            check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(nzindex)) {
        if (is.null(nzdata)) {
            nzdata <- logical(0)  # vector()
        } else if (!(is.vector(nzdata) && length(nzdata) == 0L)) {
            stop(wmsg("'nzdata' must be NULL or a vector of length 0 ",
                      "when 'nzindex' is NULL"))
        }
        nzindex <- matrix(integer(0), ncol=length(dim))
    } else {
        if (!is.matrix(nzindex))
            stop(wmsg("'nzindex' must be a matrix"))
        if (storage.mode(nzindex) == "double")
            storage.mode(nzindex) <- "integer"
        if (!is.null(dimnames(nzindex)))
            dimnames(nzindex) <- NULL
        nzdata <- .normarg_nzdata(nzdata, nrow(nzindex))
    }
    dimnames <- normarg_dimnames(dimnames, dim)
    new2("SparseArraySeed", dim=dim, nzindex=nzindex, nzdata=nzdata,
                            dimnames=dimnames,
                            check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object that supports 'type()' and subsetting
### by an M-index subscript.
### Return a SparseArraySeed object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_not_zero <- x != zero
    nzindex <- which(is_not_zero | is.na(is_not_zero), arr.ind=TRUE)  # M-index
    SparseArraySeed(x_dim, nzindex, x[nzindex], dimnames(x), check=FALSE)
}

### 'sas' must be a SparseArraySeed object.
### Return an ordinary array.
sparse2dense <- function(sas)
{
    if (!is(sas, "SparseArraySeed"))
        stop(wmsg("'sas' must be a SparseArraySeed object"))
    sas_nzdata <- nzdata(sas)
    zero <- vector(typeof(sas_nzdata), length=1L)
    ans <- array(zero, dim=dim(sas))
    ans[nzindex(sas)] <- sas_nzdata
    set_dimnames(ans, dimnames(sas))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The is_sparse() and extract_sparse_array() generics
###

### is_sparse() detects **structural** sparsity which is a global qualitative
### property of array-like object 'x' rather than a quantitative one.
### In other words it doesn't look at the data in 'x' to decide whether 'x'
### should be considered sparse or not. Said otherwise, it is NOT about
### quantitative sparsity measured by sparsity().
### IMPORTANT: Seeds for which is_sparse() returns TRUE **must** support
### extract_sparse_array().
setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

### This is the workhorse behind read_sparse_block().
### Similar to extract_array() except that:
###   (1) The extracted array data must be returned in a SparseArraySeed
###       object. Methods should always operate on the sparse representation
###       of the data and never "expand" it, that is, never turn it into a
###       dense representation for example by doing something like
###       'dense2sparse(extract_array(x, index))'. This would defeat the
###       purpose of read_sparse_block().
###   (2) It should be called only on an array-like object 'x' for which
###       'is_sparse(x)' is TRUE.
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
                      "must be an array-like object"))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse(), extract_sparse_array(), and extract_array() methods for
### SparseArraySeed objects
###

setMethod("is_sparse", "SparseArraySeed", function(x) TRUE)

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
    x_nzindex <- x@nzindex
    for (along in seq_along(ans_dim)) {
        i <- index[[along]]
        if (is.null(i))
            next
        x_nzindex[ , along] <- match(x_nzindex[ , along], i)
    }
    keep_idx <- which(!rowAnyNAs(x_nzindex))
    ans_nzindex <- x_nzindex[keep_idx, , drop=FALSE]
    ans_nzdata <- x@nzdata[keep_idx]
    SparseArraySeed(ans_dim, ans_nzindex, ans_nzdata, check=FALSE)
}
setMethod("extract_sparse_array", "SparseArraySeed",
    .extract_sparse_array_from_SparseArraySeed
)

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
    ## SparseArraySeed form (i.e. without converting first to dense then
    ## back to sparse in the process).
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
### Coercion to/from SparseArraySeed
###

### S3/S4 combo for as.array.SparseArraySeed
as.array.SparseArraySeed <- function(x, ...) sparse2dense(x)
setMethod("as.array", "SparseArraySeed", as.array.SparseArraySeed)

### S3/S4 combo for as.matrix.SparseArraySeed
.from_SparseArraySeed_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    sparse2dense(x)
}
as.matrix.SparseArraySeed <-
    function(x, ...) .from_SparseArraySeed_to_matrix(x, ...)
setMethod("as.matrix", "SparseArraySeed", .from_SparseArraySeed_to_matrix)

setAs("ANY", "SparseArraySeed", function(from) dense2sparse(from))

### Going back and forth between SparseArraySeed and dg[C|R]Matrix or
### lg[C|R]Matrix objects from the Matrix package:

.from_SparseArraySeed_to_CsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "or lgCMatrix must have exactly 2 dimensions"))
    i <- from@nzindex[ , 1L]
    j <- from@nzindex[ , 2L]
    CsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

.from_SparseArraySeed_to_RsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgRMatrix ",
                  "or lgRMatrix must have exactly 2 dimensions"))
    i <- from@nzindex[ , 1L]
    j <- from@nzindex[ , 2L]
    RsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

setAs("SparseArraySeed", "CsparseMatrix",
    .from_SparseArraySeed_to_CsparseMatrix
)
setAs("SparseArraySeed", "RsparseMatrix",
    .from_SparseArraySeed_to_RsparseMatrix
)
setAs("SparseArraySeed", "sparseMatrix",
    .from_SparseArraySeed_to_CsparseMatrix
)
setAs("SparseArraySeed", "dgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "dgCMatrix")
)
setAs("SparseArraySeed", "dgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "dgRMatrix")
)
setAs("SparseArraySeed", "lgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "lgCMatrix")
)
### Will fail if 'as(from, "RsparseMatrix")' returns a dgRMatrix object
### because the Matrix package doesn't support coercion from dgRMatrix
### to lgRMatrix at the moment:
###   > as(SparseArraySeed(4:3, rbind(c(4, 3)), -2), "lgRMatrix")
###   Error in as(as(from, "RsparseMatrix"), "lgRMatrix") :
###     no method or default for coercing “dgRMatrix” to “lgRMatrix”
setAs("SparseArraySeed", "lgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "lgRMatrix")
)

.make_SparseArraySeed_from_dgCMatrix_or_lgCMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    ans_nzindex <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    SparseArraySeed(dim(from), ans_nzindex, from@x, ans_dimnames, check=FALSE)
}

.make_SparseArraySeed_from_dgRMatrix_or_lgRMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- rep.int(seq_len(nrow(from)), diff(from@p))
    j <- from@j + 1L
    ans_nzindex <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    SparseArraySeed(dim(from), ans_nzindex, from@x, ans_dimnames, check=FALSE)
}

setAs("dgCMatrix", "SparseArraySeed",
    function(from) .make_SparseArraySeed_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("dgRMatrix", "SparseArraySeed",
    function(from) .make_SparseArraySeed_from_dgRMatrix_or_lgRMatrix(from)
)
setAs("lgCMatrix", "SparseArraySeed",
    function(from) .make_SparseArraySeed_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("lgRMatrix", "SparseArraySeed",
    function(from) .make_SparseArraySeed_from_dgRMatrix_or_lgRMatrix(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() and extract_sparse_array() methods for dg[C|R]Matrix and
### lg[C|R]Matrix objects from the Matrix package
###
### TODO: Support more sparseMatrix derivatives (e.g. dgTMatrix, dgRMatrix)
### as the need arises.
###

setMethod("is_sparse", "dgCMatrix", function(x) TRUE)
setMethod("is_sparse", "dgRMatrix", function(x) TRUE)
setMethod("is_sparse", "lgCMatrix", function(x) TRUE)
setMethod("is_sparse", "lgRMatrix", function(x) TRUE)

.extract_sparse_array_from_dgCMatrix_or_lgCMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgCMatrix or lgCMatrix object
    .make_SparseArraySeed_from_dgCMatrix_or_lgCMatrix(sm, use.dimnames=FALSE)
}

.extract_sparse_array_from_dgRMatrix_or_lgRMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgRMatrix or lgRMatrix object
    .make_SparseArraySeed_from_dgRMatrix_or_lgRMatrix(sm, use.dimnames=FALSE)
}

setMethod("extract_sparse_array", "dgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "dgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)
setMethod("extract_sparse_array", "lgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "lgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Extend base::aperm() by allowing dropping and/or adding ineffective
### dimensions. See aperm2.R
###

.aperm.SparseArraySeed <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- normarg_perm(perm, a_dim)
    msg <- validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    new_dim <- a_dim[perm]
    new_dim[is.na(perm)] <- 1L
    new_nzindex <- a@nzindex[ , perm, drop=FALSE]
    new_nzindex[ , is.na(perm)] <- 1L
    new_dimnames <- a@dimnames[perm]
    BiocGenerics:::replaceSlots(a, dim=new_dim,
                                   nzindex=new_nzindex,
                                   dimnames=new_dimnames,
                                   check=FALSE)
}

### S3/S4 combo for aperm.SparseArraySeed
aperm.SparseArraySeed <-
    function(a, perm, ...) .aperm.SparseArraySeed(a, perm, ...)
setMethod("aperm", "SparseArraySeed", aperm.SparseArraySeed)

