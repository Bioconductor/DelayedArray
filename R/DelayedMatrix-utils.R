### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() / colsum()
###

.compute_rowsum_for_block <- function(x, grid, i, j, x_is_sparse,
                                      group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    if (x_is_sparse) {
        block <- read_sparse_block(x, viewport)
        ## Go from SparseArraySeed to dgCMatrix.
        block <- as(block, "dgCMatrix")
    } else {
        block <- read_block(x, viewport)
    }
    group2 <- extractROWS(group, ranges(viewport)[1L])
    rowsum(block, group2, reorder=FALSE, na.rm=na.rm)
}
.compute_colsum_for_block <- function(x, grid, i, j, x_is_sparse,
                                      group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    if (x_is_sparse) {
        block <- read_sparse_block(x, viewport)
        ## Go from SparseArraySeed to dgCMatrix.
        block <- as(block, "dgCMatrix")
    } else {
        block <- read_block(x, viewport)
    }
    group2 <- extractROWS(group, ranges(viewport)[2L])
    colsum(block, group2, reorder=FALSE, na.rm=na.rm)
}

.compute_rowsum_for_grid_col <- function(x, grid, j, x_is_sparse,
                                         group, ugroup, na.rm=FALSE,
                                         verbose=FALSE)
{
    grid_nrow <- nrow(grid)
    grid_ncol <- ncol(grid)
    ans <- matrix(0L, nrow=length(ugroup), ncol=ncol(grid[[1L, j]]))
    ## Inner loop on the grid rows. Sequential.
    for (i in seq_len(grid_nrow)) {
        if (verbose)
            message("Processing block [[", i, "/", grid_nrow, ", ",
                                           j, "/", grid_ncol, "]] ... ",
                    appendLF=FALSE)
        block_ans <- .compute_rowsum_for_block(x, grid, i, j, x_is_sparse,
                                               group, na.rm=na.rm)
        m <- match(rownames(block_ans), ugroup)
        ans[m, ] <- ans[m, ] + block_ans
        if (verbose)
            message("OK")
    }
    ans
}
.compute_colsum_for_grid_row <- function(x, grid, i, x_is_sparse,
                                         group, ugroup, na.rm=FALSE,
                                         verbose=FALSE)
{
    grid_nrow <- nrow(grid)
    grid_ncol <- ncol(grid)
    ans <- matrix(0L, nrow=nrow(grid[[i, 1L]]), ncol=length(ugroup))
    ## Inner loop on the grid cols. Sequential.
    for (j in seq_len(grid_ncol)) {
        if (verbose)
            message("Processing block [[", i, "/", grid_nrow, ", ",
                                           j, "/", grid_ncol, "]] ... ",
                    appendLF=FALSE)
        block_ans <- .compute_colsum_for_block(x, grid, i, j, x_is_sparse,
                                               group, na.rm=na.rm)
        m <- match(colnames(block_ans), ugroup)
        ans[ , m] <- ans[ , m] + block_ans
        if (verbose)
            message("OK")
    }
    ans
}

.BLOCK_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- as.character(compute_ugroup(group, nrow(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)
    x_is_sparse <- is_sparse(x)

    ## Outer loop on the grid columns. Parallelized.
    block_results <- bplapply2(seq_len(ncol(grid)),
        function(j) {
            .compute_rowsum_for_grid_col(x, grid, j, x_is_sparse,
                                         group, ugroup, na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(cbind, block_results)
    dimnames(ans) <- list(ugroup, colnames(x))
    ans
}
.BLOCK_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- as.character(compute_ugroup(group, ncol(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)
    x_is_sparse <- is_sparse(x)

    ## Outer loop on the grid rows. Parallelized.
    block_results <- bplapply2(seq_len(nrow(grid)),
        function(i) {
            .compute_colsum_for_grid_row(x, grid, i, x_is_sparse,
                                         group, ugroup, na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(rbind, block_results)
    dimnames(ans) <- list(rownames(x), ugroup)
    ans
}

### S3/S4 combo for rowsum.DelayedMatrix
rowsum.DelayedMatrix <- function(x, group, reorder=TRUE, ...)
    .BLOCK_rowsum(x, group, reorder=reorder, ...)
setMethod("rowsum", "DelayedMatrix", .BLOCK_rowsum)

setMethod("colsum", "DelayedMatrix", .BLOCK_colsum)
