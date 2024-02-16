### =========================================================================
### rowsum() and colsum() methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###
### These methods are block processed.
###


### x: a matrix-like object (typically a DelayedMatrix).
### Walks on the grid defined on matrix-like object 'x'.
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of cbind'ing
### several pristine DelayedMatrix objects together (delayed cbind()).
BLOCK_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE,
                         grid=NULL, as.sparse=NA,
                         BACKEND=getAutoRealizationBackend(),
                         BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    stopifnot(length(dim(x)) == 2L)  # matrix-like object
    ugroup <- as.character(S4Arrays:::compute_ugroup(group, nrow(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    INIT <- function(j, grid, ugroup, x_colnames) {
        vp <- grid[[1L, j]]
        dn <- list(ugroup, extractROWS(x_colnames, ranges(vp)[2L]))
        matrix(0L, nrow=length(ugroup), ncol=ncol(vp), dimnames=dn)
    }
    INIT_MoreArgs <- list(ugroup=ugroup, x_colnames=colnames(x))

    FUN <- function(init, block, group, ugroup, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        group2 <- extractROWS(group, ranges(vp)[1L])
        block_ans <- rowsum(block, group2, reorder=FALSE, na.rm=na.rm)
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        m <- match(rownames(block_ans), ugroup)
        init[m, ] <- init[m, ] + block_ans
        init
    }
    FUN_MoreArgs <- list(group=group, ugroup=ugroup, na.rm=na.rm)

    ## No-op if 'BACKEND' is NULL.
    FINAL <- function(init, BACKEND) realize(init, BACKEND=BACKEND)
    FINAL_MoreArgs <- list(BACKEND=BACKEND)

    strip_results <- vstrip_apply(x, INIT, INIT_MoreArgs,
                                     FUN, FUN_MoreArgs,
                                     FINAL, FINAL_MoreArgs,
                                     grid=grid, as.sparse=as.sparse,
                                     BPPARAM=BPPARAM, verbose=verbose)
    do.call(cbind, strip_results)
}

### x: a matrix-like object (typically a DelayedMatrix).
### Walks on the grid defined on matrix-like object 'x'.
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of rbind'ing
### several pristine DelayedMatrix objects together (delayed rbind()).
BLOCK_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE,
                         grid=NULL, as.sparse=NA,
                         BACKEND=getAutoRealizationBackend(),
                         BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    stopifnot(length(dim(x)) == 2L)  # matrix-like object
    ugroup <- as.character(S4Arrays:::compute_ugroup(group, ncol(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    INIT <- function(i, grid, ugroup, x_rownames) {
        vp <- grid[[i, 1L]]
        dn <- list(extractROWS(x_rownames, ranges(vp)[1L]), ugroup)
        matrix(0L, nrow=nrow(vp), ncol=length(ugroup), dimnames=dn)
    }
    INIT_MoreArgs <- list(ugroup=ugroup, x_rownames=rownames(x))

    FUN <- function(init, block, group, ugroup, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        group2 <- extractROWS(group, ranges(vp)[2L])
        block_ans <- colsum(block, group2, reorder=FALSE, na.rm=na.rm)
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        m <- match(colnames(block_ans), ugroup)
        init[ , m] <- init[ , m] + block_ans
        init
    }
    FUN_MoreArgs <- list(group=group, ugroup=ugroup, na.rm=na.rm)

    ## No-op if 'BACKEND' is NULL.
    FINAL <- function(init, BACKEND) realize(init, BACKEND=BACKEND)
    FINAL_MoreArgs <- list(BACKEND=BACKEND)

    strip_results <- hstrip_apply(x, INIT, INIT_MoreArgs,
                                     FUN, FUN_MoreArgs,
                                     FINAL, FINAL_MoreArgs,
                                     grid=grid, as.sparse=as.sparse,
                                     BPPARAM=BPPARAM, verbose=verbose)
    do.call(rbind, strip_results)
}

### S3/S4 combo for rowsum.DelayedMatrix
rowsum.DelayedMatrix <- function(x, group, reorder=TRUE, ...)
    BLOCK_rowsum(x, group, reorder=reorder, ...)
setMethod("rowsum", "DelayedMatrix", BLOCK_rowsum)

setMethod("colsum", "DelayedMatrix", BLOCK_colsum)

