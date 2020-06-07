### =========================================================================
### Row and column summarization methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Workhorses behind block-processed **row** summarization functions
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

### Workhorse behind block-processed rowSums(), rowMaxs(), and rowMins().
### Could also be used for other row summarization functions (e.g. rowProds(),
### rowAnys(), rowAlls(), rowAnyNAs(), etc...). An important requirement is
### that the function must satisfy the following property: for any ordinary
### matrix 'm', '.Generic(m)' returns an atomic vector of length 'nrow(m)'.
.BLOCK_row_summary <- function(.Generic, x, na.rm=FALSE, grid=NULL)
{
    GENERIC <- match.fun(.Generic)
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)

    # TODO: either make robust to parallelization,
    # or put it in the contract that the output of
    # extract_array must support the specified GENERIC.
    coerced_GENERIC <- function(x, ...) GENERIC(as.array(x), ...)

    block_results <- blockApply(x, coerced_GENERIC, na.rm=na.rm, grid=grid)
    GENERIC(matrix(unlist(block_results, use.names=FALSE), nrow=nrow(x)))
}

### Block-processed rowMeans().
.BLOCK_rowMeans <- function(x, na.rm=FALSE, grid=NULL)
{
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
    block_results <- blockApply(x, FUN, na.rm=na.rm, grid=grid)
    combined_results <- do.call(rbind, block_results)
    row_sums <- rowSums(matrix(combined_results[ , "sums"], nrow=nrow(x)))
    row_nvals <- rowSums(matrix(combined_results[ , "nvals"], nrow=nrow(x)))
    setNames(row_sums / row_nvals, rownames(x))
}

### Block-processed rowRanges().
.BLOCK_rowRanges <- function(x, na.rm=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)

    # TODO: either make robust to parallelization,
    # or put it in the contract that the output of
    # extract_array must support rowRanges.
    coerced_rowRanges <- function(x, ...) rowRanges(as.array(x), ...)

    block_results <- blockApply(x, coerced_rowRanges, na.rm=na.rm, grid=grid)
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

.rowSums_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowSums")
    ans <- .BLOCK_row_summary(rowSums, x, na.rm=na.rm)
    setNames(ans, rownames(x))
}
setMethod("rowSums", "DelayedMatrix", .rowSums_DelayedMatrix)

.colSums_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colSums")
    ans <- .BLOCK_row_summary(rowSums, t(x), na.rm=na.rm)
    setNames(ans, colnames(x))
}
setMethod("colSums", "DelayedMatrix", .colSums_DelayedMatrix)

### row/colMeans()

.rowMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    ans <- .BLOCK_rowMeans(x, na.rm=na.rm)
    setNames(ans, rownames(x))
}
setMethod("rowMeans", "DelayedMatrix", .rowMeans_DelayedMatrix)

.colMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    ans <- .BLOCK_rowMeans(t(x), na.rm=na.rm)
    setNames(ans, colnames(x))
}
setMethod("colMeans", "DelayedMatrix", .colMeans_DelayedMatrix)


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

setGeneric("rowMaxs", signature="x")
setGeneric("colMaxs", signature="x")

.rowMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                   na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowMaxs")
    .BLOCK_row_summary(rowMaxs, x, na.rm=na.rm)
}
setMethod("rowMaxs", "DelayedMatrix", .rowMaxs_DelayedMatrix)

.colMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                   na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colMaxs")
    .BLOCK_row_summary(rowMaxs, t(x), na.rm=na.rm)
}
setMethod("colMaxs", "DelayedMatrix", .colMaxs_DelayedMatrix)

### row/colMins()

setGeneric("rowMins", signature="x")
setGeneric("colMins", signature="x")

.rowMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                   na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowMins")
    .BLOCK_row_summary(rowMins, x, na.rm=na.rm)
}
setMethod("rowMins", "DelayedMatrix", .rowMins_DelayedMatrix)

.colMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                   na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colMins")
    .BLOCK_row_summary(rowMins, t(x), na.rm=na.rm)
}
setMethod("colMins", "DelayedMatrix", .colMins_DelayedMatrix)

### row/colRanges()

.rowRanges.useAsDefault <- function(x, ...) matrixStats::rowRanges(x, ...)
setGeneric("rowRanges", signature="x",
    function(x, ...) standardGeneric("rowRanges"),
    useAsDefault=.rowRanges.useAsDefault
)
setGeneric("colRanges", signature="x")

.rowRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                     na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "rowRanges")
    .BLOCK_rowRanges(x, na.rm=na.rm)
}
setMethod("rowRanges", "DelayedMatrix", .rowRanges_DelayedMatrix)

.colRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL,
                                     na.rm=FALSE, dim.=dim(x))
{
    .check_rows_cols(rows, cols, "colRanges")
    .BLOCK_rowRanges(t(x), na.rm=na.rm)
}
setMethod("colRanges", "DelayedMatrix", .colRanges_DelayedMatrix)

### TODO: Add more row/column summarization generics/methods.

