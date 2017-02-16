### =========================================================================
### Statistical/summarization methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.normarg_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop("\"", method, "\" method for DelayedMatrix objects ",
             "does not support the 'dims' argument yet")
}

.DelayedMatrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colSums(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(m) rowSums(m, na.rm=na.rm)
    COMBINE <- function(i, m, init, reduced) { init + reduced }
    init <- numeric(nrow(x))
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm,
                                   if_empty=numeric(0))
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "DelayedMatrix", .DelayedMatrix_block_rowSums)
setMethod("colSums", "DelayedMatrix", .DelayedMatrix_block_colSums)

.DelayedMatrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "DelayedMatrix") && x@is_transposed)
        return(.DelayedMatrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(m) {
        m_sums <- rowSums(m, na.rm=na.rm)
        m_nvals <- ncol(m)
        if (na.rm)
            m_nvals <- m_nvals - rowSums(is.na(m))
        cbind(m_sums, m_nvals)
    }
    COMBINE <- function(i, m, init, reduced) { init + reduced }
    init <- cbind(
        numeric(nrow(x)),  # sums
        numeric(nrow(x))   # nvals
    )
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans[ , 1L] / ans[ , 2L], rownames(x))
}

.DelayedMatrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm,
                                    if_empty=numeric(0))
    unlist(colmeans_list, recursive=FALSE)
}

setMethod("rowMeans", "DelayedMatrix", .DelayedMatrix_block_rowMeans)
setMethod("colMeans", "DelayedMatrix", .DelayedMatrix_block_colMeans)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Row/column summarization from the matrixStats package
###
### row/colMaxs(), row/colMins(), row/colRanges(), row/colProds(),
### row/colAnys(), row/colAlls(), row/colMedians()
###

.DelayedMatrix_block_rowMaxs <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colMaxs(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))
    REDUCE <- function(m) rowMaxs(m, na.rm=na.rm)
    COMBINE <- function(i, m, init, reduced) pmax(init, reduced)
    init <- rep.int(-Inf, nrow(x))
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colMaxs <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMaxs(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))
    colmaxs_list <- colblock_APPLY(x, colMaxs, na.rm=na.rm, if_empty=-Inf)
    unlist(colmaxs_list, recursive=FALSE)
}

setGeneric("rowMaxs", signature="x")
setGeneric("colMaxs", signature="x")

setMethod("rowMaxs", "DelayedMatrix", .DelayedMatrix_block_rowMaxs)
setMethod("colMaxs", "DelayedMatrix", .DelayedMatrix_block_colMaxs)


### TODO: Add more row/column summarization generics/methods.

