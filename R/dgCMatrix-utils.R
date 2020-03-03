### =========================================================================
### Some utilities to operate on dgCMatrix objects
### -------------------------------------------------------------------------
###
### TODO: Stuff in this file should go somewhere else e.g. a better place
### would be the MatrixGenerics and sparseMatrixStats packages.
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

.dgCMatrix_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    ugroup <- compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    ans <- .Call2("C_dgCMatrix_rowsum", x, group, length(ugroup), na.rm,
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
    .dgCMatrix_rowsum(x, group, reorder=reorder, ...)
setMethod("rowsum", "dgCMatrix", rowsum.dgCMatrix)

