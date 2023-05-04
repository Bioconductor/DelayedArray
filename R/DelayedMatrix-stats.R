### =========================================================================
### Row and column summarization methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .BLOCK_row_summary()
###
### Obsolete! Use .process_grid_by_row() and .process_grid_by_col() instead.
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

### Workhorse behind rowSums(), rowMaxs(), and rowlMins() methods for
### DelayedMatrix objects.
### Could also be used for other row summarization methods e.g.
### rowProds(), rowAnys(), rowAlls(), rowAnyNAs(), etc... An important
### requirement is that '.Generic' must satisfy the following property:
### for any matrix 'm', 'match.fun(.Generic)(m)' must return an
### atomic vector of length 'nrow(m)'.
.BLOCK_row_summary <- function(.Generic, x, na.rm=FALSE, grid=NULL)
{
    GENERIC <- match.fun(.Generic)
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)
    block_results <- blockApply(x, GENERIC, na.rm=na.rm,
                                grid=grid, as.sparse=FALSE)
    GENERIC(matrix(unlist(block_results, use.names=FALSE), nrow=nrow(x)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .process_grid_by_row() and .process_grid_by_col()
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

.BLOCK_rowSums <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    INIT <- function(grid, i) vector("numeric", length=nrow(grid[[i, 1L]]))
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on rowSums() method for matrix, dgCMatrix or lgCMatrix.
        init + rowSums(block, na.rm=na.rm)
    }
    FUN_args <- list(na.rm=na.rm)

    row_summaries <-
        .process_grid_by_row(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(row_summaries, recursive=FALSE, use.names=FALSE)
}

.BLOCK_colSums <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    INIT <- function(grid, j) vector("numeric", length=ncol(grid[[1L, j]]))
    INIT_args <- list()

    FUN <- function(init, block, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        ## Dispatch on colSums() method for matrix, dgCMatrix or lgCMatrix.
        init + colSums(block, na.rm=na.rm)
    }
    FUN_args <- list(na.rm=na.rm)

    col_summaries <-
        .process_grid_by_col(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(col_summaries, recursive=FALSE, use.names=FALSE)
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
    #ans <- .BLOCK_row_summary("rowSums", x, na.rm=na.rm)
    ans <- .BLOCK_rowSums(x, na.rm=na.rm)
    setNames(ans, rownames(x))
}
setMethod("rowSums", "DelayedMatrix", .rowSums_DelayedMatrix)

### base::colSums() has a 'dims' argument. We do NOT support it.
.colSums_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colSums")
    #ans <- .BLOCK_row_summary("rowSums", t(x), na.rm=na.rm)
    ans <- .BLOCK_colSums(x, na.rm=na.rm)
    setNames(ans, colnames(x))
}
setMethod("colSums", "DelayedMatrix", .colSums_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMeans()
###
### The methods for ordinary matrices are defined in the base package.
### They propagate the rownames or colnames.
###

### TODO: Use .process_grid_by_row() instead of blockApply() to walk on 'x'.
### See .BLOCK_rowSums() above for an example.
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
    block_results <- blockApply(x, FUN, na.rm=na.rm,
                                grid=grid, as.sparse=FALSE)
    combined_results <- do.call(rbind, block_results)
    row_sums <- rowSums(matrix(combined_results[ , "sums"], nrow=nrow(x)))
    row_nvals <- rowSums(matrix(combined_results[ , "nvals"], nrow=nrow(x)))
    setNames(row_sums / row_nvals, rownames(x))
}

### base::rowMeans() has a 'dims' argument. We do NOT support it.
.rowMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    ans <- .BLOCK_rowMeans(x, na.rm=na.rm)
    setNames(ans, rownames(x))
}
setMethod("rowMeans", "DelayedMatrix", .rowMeans_DelayedMatrix)

### base::colMeans() has a 'dims' argument. We do NOT support it.
.colMeans_DelayedMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    ans <- .BLOCK_rowMeans(t(x), na.rm=na.rm)
    setNames(ans, colnames(x))
}
setMethod("colMeans", "DelayedMatrix", .colMeans_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMins()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

.BLOCK_rowMins <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
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
            block_rowmins <- MatrixGenerics::rowMins(block, na.rm=na.rm)
        }
        if (is.null(init))
            return(block_rowmins)
        pmin(init, block_rowmins)
    }
    FUN_args <- list(na.rm=na.rm)

    row_summaries <-
        .process_grid_by_row(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(row_summaries, recursive=FALSE, use.names=FALSE)
}

.BLOCK_colMins <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
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
            block_colmins <- MatrixGenerics::colMins(block, na.rm=na.rm)
        }
        if (is.null(init))
            return(block_colmins)
        pmin(init, block_colmins)
    }
    FUN_args <- list(na.rm=na.rm)

    col_summaries <-
        .process_grid_by_col(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(col_summaries, recursive=FALSE, use.names=FALSE)
}

.check_rows_cols <- function(rows, cols, method)
{
    if (!(is.null(rows) && is.null(cols)))
        stop(wmsg("\"", method, "\" method for DelayedMatrix objects ",
                  "does not support arguments 'rows' and 'cols'"))
}

### MatrixGenerics::rowMins() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "rowMins")
    #.BLOCK_row_summary("rowMins", x, na.rm=na.rm)
    .BLOCK_rowMins(x, na.rm=na.rm)
}
setMethod("rowMins", "DelayedMatrix", .rowMins_DelayedMatrix)

### MatrixGenerics::colMins() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colMins_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "colMins")
    #.BLOCK_row_summary("rowMins", t(x), na.rm=na.rm)
    .BLOCK_colMins(x, na.rm=na.rm)
}
setMethod("colMins", "DelayedMatrix", .colMins_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colMaxs()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

.BLOCK_rowMaxs <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
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
            block_rowmaxs <- MatrixGenerics::rowMaxs(block, na.rm=na.rm)
        }
        if (is.null(init))
            return(block_rowmaxs)
        pmax(init, block_rowmaxs)
    }
    FUN_args <- list(na.rm=na.rm)

    row_summaries <-
        .process_grid_by_row(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(row_summaries, recursive=FALSE, use.names=FALSE)
}

.BLOCK_colMaxs <-
    function(x, na.rm=FALSE,
             grid=NULL, as.sparse=NA, BPPARAM=getAutoBPPARAM(),
             verbose=get_verbose_block_processing())
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
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
            block_colmaxs <- MatrixGenerics::colMaxs(block, na.rm=na.rm)
        }
        if (is.null(init))
            return(block_colmaxs)
        pmax(init, block_colmaxs)
    }
    FUN_args <- list(na.rm=na.rm)

    col_summaries <-
        .process_grid_by_col(INIT, INIT_args, FUN, FUN_args, x,
                             grid=grid, as.sparse=as.sparse, BPPARAM=BPPARAM,
                             verbose=verbose)
    unlist(col_summaries, recursive=FALSE, use.names=FALSE)
}

### MatrixGenerics::rowMaxs() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "rowMaxs")
    #.BLOCK_row_summary("rowMaxs", x, na.rm=na.rm)
    .BLOCK_rowMaxs(x, na.rm=na.rm)
}
setMethod("rowMaxs", "DelayedMatrix", .rowMaxs_DelayedMatrix)

### MatrixGenerics::colMaxs() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colMaxs_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "colMaxs")
    #.BLOCK_row_summary("rowMaxs", t(x), na.rm=na.rm)
    .BLOCK_colMaxs(x, na.rm=na.rm)
}
setMethod("colMaxs", "DelayedMatrix", .colMaxs_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colRanges()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

### TODO: Use .process_grid_by_row() instead of blockApply() to walk on 'x'.
### See .BLOCK_rowMins() and .BLOCK_rowMaxs() above for examples.
.BLOCK_rowRanges <- function(x, na.rm=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(na.rm))
        stop("'na.rm' must be TRUE or FALSE")
    .get_ans_type(x, must.be.numeric=TRUE)  # only to check input type (we
                                            # ignore returned ans type)
    block_results <- blockApply(x, MatrixGenerics::rowRanges, na.rm=na.rm,
                                grid=grid, as.sparse=FALSE)
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

### MatrixGenerics::rowRanges() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.rowRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "rowRanges")
    .BLOCK_rowRanges(x, na.rm=na.rm)
}
setMethod("rowRanges", "DelayedMatrix", .rowRanges_DelayedMatrix)

### MatrixGenerics::colRanges() has the 'rows' and 'cols' arguments.
### We do NOT support them.
.colRanges_DelayedMatrix <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    .check_rows_cols(rows, cols, "colRanges")
    .BLOCK_rowRanges(t(x), na.rm=na.rm)
}
setMethod("colRanges", "DelayedMatrix", .colRanges_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### row/colVars()
###
### The methods for ordinary matrices are defined in the matrixStats package.
### They do NOT propagate the rownames or colnames.
###

### TODO


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODO: Add more row/column summarization generics/methods
###

