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

.normarg_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop("\"", method, "\" method for DelayedMatrix objects ",
             "does not support the 'dims' argument yet")
}

### row/colSums()

.DelayedMatrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colSums(t(x), na.rm=na.rm, dims=dims))

    .get_ans_type(x)  # check input type
    APPLY <- function(m) rowSums(m, na.rm=na.rm)
    COMBINE <- function(b, m, init, reduced) { init + reduced }
    init <- numeric(nrow(x))
    ans <- colblock_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    .get_ans_type(x)  # check input type
    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm)
    if (length(colsums_list) == 0L)
        return(numeric(ncol(x)))
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "DelayedMatrix", .DelayedMatrix_block_rowSums)
setMethod("colSums", "DelayedMatrix", .DelayedMatrix_block_colSums)

### row/colMeans()

.DelayedMatrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "DelayedMatrix") && x@is_transposed)
        return(.DelayedMatrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

    .get_ans_type(x)  # check input type
    APPLY <- function(m) {
        m_sums <- rowSums(m, na.rm=na.rm)
        m_nvals <- ncol(m)
        if (na.rm)
            m_nvals <- m_nvals - rowSums(is.na(m))
        cbind(m_sums, m_nvals)
    }
    COMBINE <- function(b, m, init, reduced) { init + reduced }
    init <- cbind(
        numeric(nrow(x)),  # sums
        numeric(nrow(x))   # nvals
    )
    ans <- colblock_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    setNames(ans[ , 1L] / ans[ , 2L], rownames(x))
}

.DelayedMatrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    .get_ans_type(x)  # check input type
    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm)
    if (length(colmeans_list) == 0L)
        return(rep.int(NaN, ncol(x)))
    unlist(colmeans_list, recursive=FALSE)
}

setMethod("rowMeans", "DelayedMatrix", .DelayedMatrix_block_rowMeans)
setMethod("colMeans", "DelayedMatrix", .DelayedMatrix_block_colMeans)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Row/column summarization from the matrixStats package
###
### row/colMaxs(), row/colMins(), row/colRanges(),
### row/colProds(), row/colAnys(), row/colAlls(), row/colMedians()
###

.fix_type <- function(x, ans_type)
{
    if (ans_type == "integer" && !is.integer(x) && all(is.finite(x)))
        storage.mode(x) <- ans_type
    x
}

### row/colMaxs()

.DelayedMatrix_block_rowMaxs <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colMaxs(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    APPLY <- function(m) rowMaxs(m, na.rm=na.rm)
    COMBINE <- function(b, m, init, reduced)
        .fix_type(pmax(init, reduced), ans_type)
    init <- .fix_type(rep.int(-Inf, nrow(x)), ans_type)
    ans <- colblock_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colMaxs <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMaxs(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    colmaxs_list <- colblock_APPLY(x, colMaxs, na.rm=na.rm)
    if (length(colmaxs_list) == 0L)
        return(.fix_type(rep.int(-Inf, ncol(x)), ans_type))
    unlist(colmaxs_list, recursive=FALSE)
}

setGeneric("rowMaxs", signature="x")
setGeneric("colMaxs", signature="x")

setMethod("rowMaxs", "DelayedMatrix", .DelayedMatrix_block_rowMaxs)
setMethod("colMaxs", "DelayedMatrix", .DelayedMatrix_block_colMaxs)

### row/colMins()

.DelayedMatrix_block_rowMins <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colMins(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    APPLY <- function(m) rowMins(m, na.rm=na.rm)
    COMBINE <- function(b, m, init, reduced)
        .fix_type(pmin(init, reduced), ans_type)
    init <- .fix_type(rep.int(Inf, nrow(x)), ans_type)
    ans <- colblock_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colMins <- function(x, rows=NULL, cols=NULL,
                                         na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMins(t(x), rows=rows, cols=cols,
                                            na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    colmins_list <- colblock_APPLY(x, colMins, na.rm=na.rm)
    if (length(colmins_list) == 0L)
        return(.fix_type(rep.int(Inf, ncol(x)), ans_type))
    unlist(colmins_list, recursive=FALSE)
}

setGeneric("rowMins", signature="x")
setGeneric("colMins", signature="x")

setMethod("rowMins", "DelayedMatrix", .DelayedMatrix_block_rowMins)
setMethod("colMins", "DelayedMatrix", .DelayedMatrix_block_colMins)

### row/colRanges()

.DelayedMatrix_block_rowRanges <- function(x, rows=NULL, cols=NULL,
                                           na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colRanges(t(x), rows=rows, cols=cols,
                                              na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    APPLY <- function(m) rowRanges(m, na.rm=na.rm)
    COMBINE <- function(b, m, init, reduced) {
        .fix_type(cbind(pmin(init[ , 1L], reduced[ , 1L]),
                        pmax(init[ , 2L], reduced[ , 2L])),
                  ans_type)
    }
    init <- .fix_type(matrix(rep(c(Inf, -Inf), each=nrow(x)), ncol=2L),
                      ans_type)
    ans <- colblock_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colRanges <- function(x, rows=NULL, cols=NULL,
                                           na.rm=FALSE, dim.=dim(x))
{
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowRanges(t(x), rows=rows, cols=cols,
                                              na.rm=na.rm, dim.=dim.))

    ans_type <- .get_ans_type(x, must.be.numeric=TRUE)
    colranges_list <- colblock_APPLY(x, colRanges, na.rm=na.rm)
    if (length(colranges_list) == 0L)
        return(.fix_type(matrix(rep(c(Inf, -Inf), each=ncol(x)), ncol=2L),
                         ans_type))
    do.call(rbind, colranges_list)
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

