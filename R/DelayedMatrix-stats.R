### =========================================================================
### Statistical/summarization methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###


### Raise an error if invalid input type. Otherwise return "integer",
### "numeric", "double", or "complex".
.get_ans_type <- function(x, must.be.numeric=FALSE)
{
    x_type <- type(x)
    ans_type <- switch(x_type,
        logical="integer",
        integer=, numeric=, double=, complex=x_type,
        stop(wmsg("operation not supported on matrices of type ", x_type)))
    if (must.be.numeric && !is.numeric(get(ans_type)(0)))
        stop(wmsg("operation not supported on matrices of type ", x_type))
    ans_type
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colSums() and row/colMeans()
###

.check_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop(wmsg("\"", method, "\" method for DelayedMatrix objects ",
                  "does not support the 'dims' argument"))
}

### row/colSums()

.DelayedMatrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowSums")
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x)  # only to check input type (we ignore returned ans type)

    block_results <- blockApply(x, rowSums, na.rm=na.rm)

    ans <- rowSums(matrix(unlist(block_results, use.names=FALSE), nrow=nrow(x)))
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colSums")
    .DelayedMatrix_block_rowSums(t(x), na.rm=na.rm, dims=dims)
}

setMethod("rowSums", "DelayedMatrix", .DelayedMatrix_block_rowSums)
setMethod("colSums", "DelayedMatrix", .DelayedMatrix_block_colSums)

### row/colMeans()

.DelayedMatrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x)  # only to check input type (we ignore returned ans type)

    FUN <- function(block, na.rm=FALSE) {
        block_sums <- rowSums(block, na.rm=na.rm)
        block_nvals <- rep.int(ncol(block), nrow(block))
        if (na.rm)
            block_nvals <- block_nvals - rowSums(is.na(block))
        cbind(sums=block_sums, nvals=block_nvals)
    }
    block_results <- blockApply(x, FUN, na.rm=na.rm)

    combined_results <- do.call(rbind, block_results)
    row_sums <- rowSums(matrix(combined_results[ , "sums"], nrow=nrow(x)))
    row_nvals <- rowSums(matrix(combined_results[ , "nvals"], nrow=nrow(x)))
    setNames(row_sums / row_nvals, rownames(x))
}

.DelayedMatrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    .DelayedMatrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims)
}

setMethod("rowMeans", "DelayedMatrix", .DelayedMatrix_block_rowMeans)
setMethod("colMeans", "DelayedMatrix", .DelayedMatrix_block_colMeans)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Row/column summarization functions from the matrixStats package
###
### row/colMaxs(), row/colMins(), row/colRanges(),
### row/colProds(), row/colAnys(), row/colAlls(), row/colMedians()
###
### All these functions have the 'rows', 'cols', and 'dim.' arguments. We
### ignore these arguments for now.
### Unlike row/colSums() and row/colMeans() from the base package, these
### functions don't propagate the rownames or colnames.
###

.check_rows_cols <- function(rows, cols, method)
{
    if (!(is.null(rows) && is.null(cols)))
        stop(wmsg("\"", method, "\" method for DelayedMatrix objects ",
                  "does not support arguments 'rows' and 'cols'"))
}

### row/colMaxs()

.DelayedMatrix_block_rowMaxs <- function(x, rows=NULL, cols=NULL,
                                            na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowMaxs")
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)

    block_results <- blockApply(x, rowMaxs, na.rm=na.rm)

    rowMaxs(matrix(unlist(block_results, use.names=FALSE), nrow=nrow(x)))
}

.DelayedMatrix_block_colMaxs <- function(x, rows=NULL, cols=NULL,
                                            na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colMaxs")
    .DelayedMatrix_block_rowMaxs(t(x), rows=cols, cols=rows,
                                       na.rm=na.rm, dim.=dim.)
}

setGeneric("rowMaxs", signature="x")
setGeneric("colMaxs", signature="x")

setMethod("rowMaxs", "DelayedMatrix", .DelayedMatrix_block_rowMaxs)
setMethod("colMaxs", "DelayedMatrix", .DelayedMatrix_block_colMaxs)

### row/colMins()

.DelayedMatrix_block_rowMins <- function(x, rows=NULL, cols=NULL,
                                            na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowMins")
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)

    block_results <- blockApply(x, rowMins, na.rm=na.rm)

    rowMins(matrix(unlist(block_results, use.names=FALSE), nrow=nrow(x)))
}

.DelayedMatrix_block_colMins <- function(x, rows=NULL, cols=NULL,
                                            na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colMins")
    .DelayedMatrix_block_rowMins(t(x), rows=cols, cols=rows,
                                       na.rm=na.rm, dim.=dim.)
}

setGeneric("rowMins", signature="x")
setGeneric("colMins", signature="x")

setMethod("rowMins", "DelayedMatrix", .DelayedMatrix_block_rowMins)
setMethod("colMins", "DelayedMatrix", .DelayedMatrix_block_colMins)

### row/colRanges()

.DelayedMatrix_block_rowRanges <- function(x, rows=NULL, cols=NULL,
                                              na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowRanges")
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)

    block_results <- blockApply(x, rowRanges, na.rm=na.rm)

    combined_results <- do.call(rbind, block_results)
    row_mins <- rowMins(matrix(combined_results[ , 1L], nrow=nrow(x)))
    row_maxs <- rowMaxs(matrix(combined_results[ , 2L], nrow=nrow(x)))
    ans <- cbind(row_mins, row_maxs, deparse.level=0)
    ## 'ans' can have unexpected dimnames because of the following bug
    ## in cbind/rbind:
    ##   https://stat.ethz.ch/pipermail/r-devel/2017-December/075288.html
    ## TODO: Remove the line below once the above bug is fixed.
    dimnames(ans) <- NULL
    ans
}

.DelayedMatrix_block_colRanges <- function(x, rows=NULL, cols=NULL,
                                              na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colRanges")
    .DelayedMatrix_block_rowRanges(t(x), rows=cols, cols=rows,
                                         na.rm=na.rm, dim.=dim.)
}

.rowRanges.useAsDefault <- function(x, ...) matrixStats::rowRanges(x, ...)
setGeneric("rowRanges", signature="x",
    function(x, ...) standardGeneric("rowRanges"),
    useAsDefault=.rowRanges.useAsDefault
)

setGeneric("colRanges", signature="x")

setMethod("rowRanges", "DelayedMatrix", .DelayedMatrix_block_rowRanges)
setMethod("colRanges", "DelayedMatrix", .DelayedMatrix_block_colRanges)

### TODO: Add more row/column summarization generics/methods.

