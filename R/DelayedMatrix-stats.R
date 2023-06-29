### =========================================================================
### Row and column summarization methods for DelayedMatrix objects
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
### .process_grid_by_row() and .process_grid_by_col()
###
### Obsolete! Use block_apply_by_grid_row() and block_apply_by_grid_col()
### instead.
###

.process_grid_by_row <-
    function(INIT, INIT_args, FUN, FUN_args, x,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be a FALSE, TRUE, or NA"))

    process_grid_row <- function(init, FUN, FUN_args, x,
                                 grid, i, as.sparse, verbose)
    {
        grid_nrow <- nrow(grid)
        grid_ncol <- ncol(grid)
        ## Inner loop on the grid columns. Sequential.
        for (j in seq_len(grid_ncol)) {
            if (verbose)
                message("Processing block [[", i, "/", grid_nrow, ", ",
                                               j, "/", grid_ncol, "]] ... ",
                        appendLF=FALSE)
            viewport <- grid[[i, j]]
            block <- read_block(x, viewport, as.sparse=as.sparse)
            init <- do.call(FUN, c(list(init, block), FUN_args))
            if (verbose)
                message("OK")
        }
        init
    }

    ## Outer loop on the grid rows. Parallelized.
    S4Arrays:::bplapply2(seq_len(nrow(grid)),
        function(i) {
            init <- do.call(INIT, c(list(grid, i), INIT_args))
            process_grid_row(init, FUN, FUN_args, x,
                             grid, i, as.sparse, verbose)
        },
        BPPARAM=BPPARAM
    )
}

.process_grid_by_col <-
    function(INIT, INIT_args, FUN, FUN_args, x,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be a FALSE, TRUE, or NA"))

    process_grid_col <- function(init, FUN, FUN_args, x,
                                 grid, j, as.sparse, verbose)
    {
        grid_nrow <- nrow(grid)
        grid_ncol <- ncol(grid)
        ## Inner loop on the grid rows. Sequential.
        for (i in seq_len(grid_nrow)) {
            if (verbose)
                message("Processing block [[", i, "/", grid_nrow, ", ",
                                               j, "/", grid_ncol, "]] ... ",
                        appendLF=FALSE)
            viewport <- grid[[i, j]]
            block <- read_block(x, viewport, as.sparse=as.sparse)
            init <- do.call(FUN, c(list(init, block), FUN_args))
            if (verbose)
                message("OK")
        }
        init
    }

    ## Outer loop on the grid columns. Parallelized.
    S4Arrays:::bplapply2(seq_len(ncol(grid)),
        function(j) {
            init <- do.call(INIT, c(list(grid, j), INIT_args))
            process_grid_col(init, FUN, FUN_args, x,
                             grid, j, as.sparse, verbose)
        },
        BPPARAM=BPPARAM
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colSums()
###
### The methods for ordinary matrices are defined in the base package.
### They propagate the rownames or colnames.
###

BLOCK_rowSums <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    INIT <- function(grid, i) numeric(nrow(grid[[i, 1L]]))
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on rowSums() method for matrix, dgCMatrix or lgCMatrix.
        init + rowSums(block, na.rm=na.rm)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- rownames(x)
    ans
}

BLOCK_colSums <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    INIT <- function(grid, j) numeric(ncol(grid[[1L, j]]))
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on colSums() method for matrix, dgCMatrix or lgCMatrix.
        init + colSums(block, na.rm=na.rm)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_col(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

.check_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop(wmsg("\"", method, "\" method for DelayedMatrix objects ",
                  "does not support the 'dims' argument"))
}

### base::rowSums() has a 'dims' argument. We do NOT support it.
.rowSums_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowSums")
    BLOCK_rowSums(x, na.rm=na.rm)
}
setMethod("rowSums", "DelayedMatrix", .rowSums_DelayedMatrix)

### base::colSums() has a 'dims' argument. We do NOT support it.
.colSums_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colSums")
    BLOCK_colSums(x, na.rm=na.rm)
}
setMethod("colSums", "DelayedMatrix", .colSums_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMeans()
###
### The methods for ordinary matrices are defined in the base package.
### They propagate the rownames or colnames.
###

BLOCK_rowMeans <- function(x, na.rm=FALSE, useNames=TRUE,
                           grid=NULL, as.sparse=NA,
                           BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    INIT <- function(grid, i) {
        n <- nrow(grid[[i, 1L]])
        data.frame(sums=numeric(n), nvals=integer(n))
    }
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on rowSums() method for matrix, dgCMatrix or lgCMatrix.
        block_sums <- rowSums(block, na.rm=na.rm)
        block_nvals <- rep.int(ncol(block), nrow(block))
        if (na.rm)
            block_nvals <- block_nvals - rowSums(is.na(block))
        init + data.frame(sums=block_sums, nvals=block_nvals)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    df <- do.call(rbind, res_list)
    ans <- df$sums / df$nvals
    if (useNames)
        names(ans) <- rownames(x)
    ans
}

BLOCK_colMeans <- function(x, na.rm=FALSE, useNames=TRUE,
                           grid=NULL, as.sparse=NA,
                           BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    INIT <- function(grid, j) {
        n <- ncol(grid[[1L, j]])
        data.frame(sums=numeric(n), nvals=integer(n))
    }
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on colSums() method for matrix, dgCMatrix or lgCMatrix.
        block_sums <- colSums(block, na.rm=na.rm)
        block_nvals <- rep.int(nrow(block), ncol(block))
        if (na.rm)
            block_nvals <- block_nvals - colSums(is.na(block))
        init + data.frame(sums=block_sums, nvals=block_nvals)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_col(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    df <- do.call(rbind, res_list)
    ans <- df$sums / df$nvals
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

### base::rowMeans() has a 'dims' argument. We do NOT support it.
.rowMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    BLOCK_rowMeans(x, na.rm=na.rm)
}
setMethod("rowMeans", "DelayedMatrix", .rowMeans_DelayedMatrix)

### base::colMeans() has a 'dims' argument. We do NOT support it.
.colMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    BLOCK_colMeans(x, na.rm=na.rm)
}
setMethod("colMeans", "DelayedMatrix", .colMeans_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMins()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

BLOCK_rowMins <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    if (ncol(x) == 0L)
        return(rep.int(Inf, nrow(x)))

    INIT <- function(grid, i) NULL
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed")) {
            ## Transposing a SparseArraySeed object is faster than transposing
            ## a CsparseMatrix derivative so we transpose **before** coercion
            ## to CsparseMatrix.
            block <- as(t(block), "CsparseMatrix")  # to dgCMatrix or lgCMatrix
            block_rowmins <- SparseArray:::colMins_dgCMatrix(block, na.rm=na.rm)
        } else {
            ## Dispatch on matrixStats::rowMins().
            block_rowmins <- MatrixGenerics::rowMins(block, na.rm=na.rm,
                                                     useNames=FALSE)
        }
        if (is.null(init))
            return(block_rowmins)
        pmin(init, block_rowmins)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- rownames(x)
    ans
}

BLOCK_colMins <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    if (nrow(x) == 0L)
        return(rep.int(Inf, ncol(x)))

    INIT <- function(grid, j) NULL
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed")) {
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
            block_colmins <- SparseArray:::colMins_dgCMatrix(block, na.rm=na.rm)
        } else {
            ## Dispatch on matrixStats::colMins().
            block_colmins <- MatrixGenerics::colMins(block, na.rm=na.rm,
                                                     useNames=FALSE)
        }
        if (is.null(init))
            return(block_colmins)
        pmin(init, block_colmins)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_col(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

.check_rows_cols <- function(rows, cols, method)
{
    if (!(is.null(rows) && is.null(cols)))
        stop(wmsg("\"", method, "\" method for DelayedMatrix objects ",
                  "does not support arguments 'rows' and 'cols'"))
}

### MatrixGenerics::rowMins() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE,
                                   useNames=TRUE)
{
    .check_rows_cols(rows, cols, "rowMins")
    BLOCK_rowMins(x, na.rm=na.rm, useNames=useNames)
}
setMethod("rowMins", "DelayedMatrix", .rowMins_DelayedMatrix)

### MatrixGenerics::colMins() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE,
                                   useNames=TRUE)
{
    .check_rows_cols(rows, cols, "colMins")
    BLOCK_colMins(x, na.rm=na.rm, useNames=useNames)
}
setMethod("colMins", "DelayedMatrix", .colMins_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMaxs()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

BLOCK_rowMaxs <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    if (ncol(x) == 0L)
        return(rep.int(-Inf, nrow(x)))

    INIT <- function(grid, i) NULL
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed")) {
            ## Transposing a SparseArraySeed object is faster than transposing
            ## a CsparseMatrix derivative so we transpose **before** coercion
            ## to CsparseMatrix.
            block <- as(t(block), "CsparseMatrix")  # to dgCMatrix or lgCMatrix
            block_rowmaxs <- SparseArray:::colMaxs_dgCMatrix(block, na.rm=na.rm)
        } else {
            ## Dispatch on matrixStats::rowMaxs().
            block_rowmaxs <- MatrixGenerics::rowMaxs(block, na.rm=na.rm,
                                                     useNames=FALSE)
        }
        if (is.null(init))
            return(block_rowmaxs)
        pmax(init, block_rowmaxs)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- rownames(x)
    ans
}

BLOCK_colMaxs <- function(x, na.rm=FALSE, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    if (nrow(x) == 0L)
        return(rep.int(-Inf, ncol(x)))

    INIT <- function(grid, j) NULL
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed")) {
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
            block_colmaxs <- SparseArray:::colMaxs_dgCMatrix(block, na.rm=na.rm)
        } else {
            ## Dispatch on matrixStats::colMaxs().
            block_colmaxs <- MatrixGenerics::colMaxs(block, na.rm=na.rm,
                                                     useNames=FALSE)
        }
        if (is.null(init))
            return(block_colmaxs)
        pmax(init, block_colmaxs)
    }
    FUN_args <- list(na.rm=na.rm)

    res_list <- block_apply_by_grid_col(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

### MatrixGenerics::rowMaxs() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE,
                                   useNames=TRUE)
{
    .check_rows_cols(rows, cols, "rowMaxs")
    BLOCK_rowMaxs(x, na.rm=na.rm, useNames=useNames)
}
setMethod("rowMaxs", "DelayedMatrix", .rowMaxs_DelayedMatrix)

### MatrixGenerics::colMaxs() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE,
                                   useNames=TRUE)
{
    .check_rows_cols(rows, cols, "colMaxs")
    BLOCK_colMaxs(x, na.rm=na.rm, useNames=useNames)
}
setMethod("colMaxs", "DelayedMatrix", .colMaxs_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colRanges()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

### TODO:
### - Use block_apply_by_grid_row() instead of blockApply() to walk on 'x'.
### - Support arguments 'useNames', 'as.sparse', 'BPPARAM', and 'verbose'.
### See BLOCK_rowMins() and BLOCK_rowMaxs() above for examples.
BLOCK_rowRanges <- function(x, na.rm=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)
    block_results <- blockApply(x, MatrixGenerics::rowRanges,
                                   na.rm=na.rm, useNames=FALSE,
                                   grid=grid, as.sparse=FALSE)
    combined_results <- do.call(rbind, block_results)
    row_mins <- rowMins(matrix(combined_results[ , 1L], nrow=nrow(x)),
                        useNames=FALSE)
    row_maxs <- rowMaxs(matrix(combined_results[ , 2L], nrow=nrow(x)),
                        useNames=FALSE)
    ans <- cbind(row_mins, row_maxs, deparse.level=0)
    ## 'ans' can have unexpected dimnames because of the following bug
    ## in cbind/rbind:
    ##   https://stat.ethz.ch/pipermail/r-devel/2017-December/075288.html
    ## TODO: Remove the line below once the above bug is fixed.
    dimnames(ans) <- NULL
    ans
}

### MatrixGenerics::rowRanges() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "rowRanges")
    BLOCK_rowRanges(x, na.rm=na.rm)
}
setMethod("rowRanges", "DelayedMatrix", .rowRanges_DelayedMatrix)

### MatrixGenerics::colRanges() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "colRanges")
    BLOCK_rowRanges(t(x), na.rm=na.rm)
}
setMethod("colRanges", "DelayedMatrix", .colRanges_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colVars()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

BLOCK_rowVars <- function(x, na.rm=FALSE, center=NULL, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    x_nrow <- nrow(x)
    if (is.null(center)) {
        center <- BLOCK_rowMeans(x, na.rm=na.rm, useNames=FALSE,
                                 grid=grid, as.sparse=as.sparse,
                                 BPPARAM=BPPARAM, verbose=verbose)
    } else {
        if (!is.numeric(center))
            stop("'center' must be NULL or a numeric vector")
        if (length(center) != x_nrow) {
            if (length(center) != 1L)
                stop("'center' must have length 1 or nrow(x)")
            center <- rep.int(center, x_nrow)
        }
    }

    INIT <- function(grid, i) {
        n <- nrow(grid[[i, 1L]])
        data.frame(sums=numeric(n), nvals=integer(n))
    }
    INIT_args <- list()

    FUN <- function(init, block, center, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        viewport <- currentViewport()
        block_center <- extractROWS(center, ranges(viewport)[1L])
        delta <- block - block_center
        ## Dispatch on rowSums() method for matrix, dgCMatrix or lgCMatrix.
        block_sums <- rowSums(delta * delta, na.rm=na.rm)
        block_nvals <- rep.int(ncol(block), nrow(block))
        if (na.rm)
            block_nvals <- block_nvals - rowSums(is.na(block))
        init + data.frame(sums=block_sums, nvals=block_nvals)
    }
    FUN_args <- list(center, na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    df <- do.call(rbind, res_list)
    ans <- df$sums / (df$nvals - 1L)
    if (useNames)
        names(ans) <- rownames(x)
    ans
}

BLOCK_colVars <- function(x, na.rm=FALSE, center=NULL, useNames=TRUE,
                          grid=NULL, as.sparse=NA,
                          BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(useNames))
        stop(wmsg("'useNames' must be TRUE or FALSE"))

    x_ncol <- ncol(x)
    if (is.null(center)) {
        center <- BLOCK_colMeans(x, na.rm=na.rm, useNames=FALSE,
                                 grid=grid, as.sparse=as.sparse,
                                 BPPARAM=BPPARAM, verbose=verbose)
    } else {
        if (!is.numeric(center))
            stop("'center' must be NULL or a numeric vector")
        if (length(center) != x_ncol) {
            if (length(center) != 1L)
                stop("'center' must have length 1 or ncol(x)")
            center <- rep.int(center, x_ncol)
        }
    }

    INIT <- function(grid, j) {
        n <- ncol(grid[[1L, j]])
        data.frame(sums=numeric(n), nvals=integer(n))
    }
    INIT_args <- list()

    FUN <- function(init, block, center, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        viewport <- currentViewport()
        block_center <- extractROWS(center, ranges(viewport)[2L])
        delta <- block - rep(block_center, each=nrow(block))
        ## Dispatch on rowSums() method for matrix, dgCMatrix or lgCMatrix.
        block_sums <- colSums(delta * delta, na.rm=na.rm)
        block_nvals <- rep.int(nrow(block), ncol(block))
        if (na.rm)
            block_nvals <- block_nvals - colSums(is.na(block))
        init + data.frame(sums=block_sums, nvals=block_nvals)
    }
    FUN_args <- list(center, na.rm=na.rm)

    res_list <- block_apply_by_grid_row(INIT, INIT_args, FUN, FUN_args, x,
                                        grid=grid, as.sparse=as.sparse,
                                        BPPARAM=BPPARAM, verbose=verbose)
    df <- do.call(rbind, res_list)
    ans <- df$sums / (df$nvals - 1L)
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

setMethod("rowVars", "DelayedMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, useNames=TRUE)
    {
        .check_rows_cols(rows, cols, "rowVars")
        BLOCK_rowVars(x, na.rm=na.rm, center=center, useNames=useNames)
    }
)

setMethod("colVars", "DelayedMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, useNames=TRUE)
    {
        .check_rows_cols(rows, cols, "colVars")
        BLOCK_colVars(x, na.rm=na.rm, center=center, useNames=useNames)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODO: Add more row/column summarization generics/methods
###

