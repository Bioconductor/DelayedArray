### =========================================================================
### Some summarization methods that operate natively on dgCMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparseMatrix2() -- NOT exported
###
### A replacement for Matrix::sparseMatrix() that is typically 50%-60% faster
### and more memory efficient. Like Matrix::sparseMatrix(), it only supports
### numeric or logical input data at the moment. If 'is.numeric(nzdata)' is
### TRUE, it returns a dgCMatrix object. If 'is.logical(nzdata)' is TRUE, it
### returns a lgCMatrix object. Any other type of input triggers an error.

### 'i', 'j', 'nzdata' must be **parallel** atomic vectors (integer vectors
### with no NAs for 'i' and 'j', and integer, double or logical vector for
### 'nzdata', possibly with NAs).
sparseMatrix2 <- function(dim, i, j, nzdata, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L,
              is.integer(i), is.integer(j))
    nzdata_type <- typeof(nzdata)
    ans_class <- switch(nzdata_type,
                        'integer'=, 'double'="dgCMatrix",
                        'logical'="lgCMatrix",
                        stop(wmsg("unsupported data type: ", nzdata_type)))
    dimnames <- normarg_dimnames(dimnames, dim)
    oo <- order(j, i)
    ans_i <- i[oo] - 1L  # dgCMatrix and lgCMatrix objects want this zero-based
    ans_p <- c(0L, cumsum(tabulate(j[oo], nbins=dim[[2L]])))
    ans_x <- nzdata[oo]
    if (is.integer(ans_x))
        ans_x <- as.double(ans_x)
    new(ans_class, Dim=dim, i=ans_i, p=ans_p, x=ans_x, Dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() and colsum() methods
###

compute_ugroup <- function(group, expected_group_len, reorder=TRUE)
{
    if (!(is.vector(group) || is.factor(group)))
        stop(wmsg("'group' must be a vector or factor"))
    if (length(group) != expected_group_len)
        stop(wmsg("incorrect length for 'group'"))
    if (!isTRUEorFALSE(reorder))
        stop(wmsg("'reorder' must be TRUE or FALSE"))
    ## Taken from base::rowsum.default().
    ugroup <- unique(group)
    if (anyNA(ugroup))
        warning(wmsg("missing values for 'group'"))
    if (reorder)
        ugroup <- sort(ugroup, na.last=TRUE, method="quick")
    ugroup
}

.rowsum_dgCMatrix <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    ugroup <- compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    ans <- .Call2("C_rowsum_dgCMatrix", x, group, length(ugroup), na.rm,
                                        PACKAGE="DelayedArray")
    dimnames(ans) <- list(as.character(ugroup), colnames(x))
    ans
}

### The base package provides rowsum() only (as an S3 generic).
setGeneric("rowsum", signature="x")

setGeneric("colsum", signature="x",
    function(x, group, reorder=TRUE, ...)
        standardGeneric("colsum")
)

setMethod("colsum", "ANY",
    function(x, group, reorder=TRUE, ...)
    {
        t(rowsum(t(x), group, reorder=reorder, ...))
    }
)

### S3/S4 combo for rowsum.dgCMatrix
rowsum.dgCMatrix <- function(x, group, reorder=TRUE, ...)
    .rowsum_dgCMatrix(x, group, reorder=reorder, ...)
setMethod("rowsum", "dgCMatrix", rowsum.dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins_dgCMatrix()
### colMaxs_dgCMatrix()
### colRanges_dgCMatrix()
### colVars_dgCMatrix()
###
### NOT exported.
###
### Don't turn these into formal S4 methods for dgCMatrix objects to avoid
### conflict with the methods defined in the sparseMatrixStats package!
### They do NOT propagate the colnames (the methods defined in matrixStats
### don't either).

colMins_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colMins_dgCMatrix", x, na.rm, PACKAGE="DelayedArray")
}

colMaxs_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colMaxs_dgCMatrix", x, na.rm, PACKAGE="DelayedArray")
}

### About 2x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colRanges_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colRanges_dgCMatrix", x, na.rm, PACKAGE="DelayedArray")
}

### About 2.5x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colVars_dgCMatrix <- function(x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colVars_dgCMatrix", x, na.rm, PACKAGE="DelayedArray")
}

