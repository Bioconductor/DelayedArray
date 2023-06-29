### =========================================================================
### blockApply() and family
### -------------------------------------------------------------------------
###


### NOT exported but used in the HDF5Array package!
get_verbose_block_processing <- function()
{
    getOption("DelayedArray.verbose.block.processing", default=FALSE)
}

### NOT exported but used in the HDF5Array package!
set_verbose_block_processing <- function(verbose)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    old_verbose <- get_verbose_block_processing()
    options(DelayedArray.verbose.block.processing=verbose)
    old_verbose
}

normarg_verbose <- function(verbose)
{
    if (!(is.logical(verbose) && length(verbose) == 1L))
        stop(wmsg("'verbose' must be FALSE, TRUE, or NA"))
    if (is.na(verbose))
        verbose <- get_verbose_block_processing()
    verbose
}

.realize_what_as_what <- function(x_is_sparse, as.sparse)
{
    if (is.na(as.sparse) || as.sparse == x_is_sparse) {
        what <- if (x_is_sparse) "sparse block" else "block"
        as_what <- ""
    } else {
        if (x_is_sparse) {
            what <- "sparse block"
            as_what <- "dense block"
        } else {
            what <- "dense block"
            as_what <- "sparse block"
        }
        as_what <- paste0(" as ", as_what)
    }
    list(what=what, as_what=as_what)
}

### For use in blockApply() and family.
verbose_read_block <- function(x, viewport, x_is_sparse, as.sparse, bid, nblock)
{
    what_as_what <- .realize_what_as_what(x_is_sparse, as.sparse)
    message("/ Reading and realizing ", what_as_what$what, " ",
            bid, "/", nblock, what_as_what$as_what, " ... ", appendLF=FALSE)
    block <- read_block(x, viewport, as.sparse=as.sparse)
    message("OK")
    block
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoBPPARAM()
###

### By default (i.e. when no argument is specified), no BiocParallel backend
### is set and evaluation is sequential.
### Beware that SnowParam() on Windows is quite inefficient for block
### processing (it introduces **a lot** of overhead) so it's better to stick
### to sequential evaluation on this platform.
### See https://github.com/Bioconductor/BiocParallel/issues/78
setAutoBPPARAM <- function(BPPARAM=NULL)
{
    if (!is.null(BPPARAM)) {
        if (!requireNamespace("BiocParallel", quietly=TRUE))
            stop(wmsg("Couldn't load the BiocParallel package. Please ",
                      "install the BiocParallel package and try again."))
        if (!is(BPPARAM, "BiocParallelParam"))
            stop(wmsg("'BPPARAM' must be a BiocParallelParam ",
                      "object from the BiocParallel package"))
    }
    S4Arrays:::set_user_option("auto.BPPARAM", BPPARAM)
}

getAutoBPPARAM <- function() S4Arrays:::get_user_option("auto.BPPARAM")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set/get grid context for the current block of a gridApply(), blockApply(),
### gridReduce(), or blockReduce() loop
###

set_grid_context <- function(effective_grid,
                             current_block_id,
                             current_viewport,
                             envir=parent.frame(1))
{
    if (!is.null(effective_grid))
        assign(".effective_grid", effective_grid, envir=envir)
    if (!is.null(current_block_id))
        assign(".current_block_id", current_block_id, envir=envir)
    if (!is.null(current_viewport))
        assign(".current_viewport", current_viewport, envir=envir)
}

.backward_compat <- function(envir, funname)
{
    if (!(is.array(envir) || is(envir, "SparseArraySeed")))
        return(envir)
    msg <- c("starting with DelayedArray 0.21.1, passing 'block' ",
             "to ", funname, "() is no longer needed and is ",
             "considered an error")
    .Defunct(msg=c("  ", wmsg(msg)))
    parent.frame(3)
}

.grid_context_not_found <- c(
    "Grid context not found for the current block. ",
    "Are we in a blockApply(), blockReduce(), gridApply(), ",
    "or gridReduce() loop?"
)

.explain_proper_use <- function(funname)
    paste0("Note that ", funname, "() can only be called from **within** ",
           "the callback functions 'FUN' and/or 'BREAKIF' passed to ",
           "blockApply() and family.")

.suggest_set_grid_context <- c(
    "If you need to be able to test/debug your callback function 'FUN' ",
    "(or 'BREAKIF') as a standalone function, set an arbitrary grid context ",
    "by calling set_grid_context() **right before** calling the callback ",
    "function."
)

effectiveGrid <- function(envir=parent.frame(2))
{
    envir <- .backward_compat(envir, "effectiveGrid")
    effective_grid <- try(get(".effective_grid", envir=envir,
                              inherits=FALSE), silent=TRUE)
    if (inherits(effective_grid, "try-error"))
        stop(wmsg(.grid_context_not_found),
             "\n\n  ",
             wmsg(.explain_proper_use("effectiveGrid")),
             "\n\n  ",
             wmsg(.suggest_set_grid_context))
    effective_grid
}

currentBlockId <- function(envir=parent.frame(2))
{
    envir <- .backward_compat(envir, "currentBlockId")
    current_block_id <- try(get(".current_block_id", envir=envir,
                                inherits=FALSE), silent=TRUE)
    if (inherits(current_block_id, "try-error"))
        stop(wmsg(.grid_context_not_found),
             "\n\n  ",
             wmsg(.explain_proper_use("currentBlockId")),
             "\n\n  ",
             wmsg(.suggest_set_grid_context))
    current_block_id
}

currentViewport <- function(envir=parent.frame(2))
{
    envir <- .backward_compat(envir, "currentViewport")
    current_viewport <- try(get(".current_viewport", envir=envir,
                                inherits=FALSE), silent=TRUE)
    if (!inherits(current_viewport, "try-error"))
        return(current_viewport)
    effective_grid <- try(effectiveGrid(envir), silent=TRUE)
    if (inherits(effective_grid, "try-error"))
        stop(wmsg(.grid_context_not_found),
             "\n\n  ",
             wmsg(.explain_proper_use("currentViewport")),
             "\n\n  ",
             wmsg(.suggest_set_grid_context))
    effective_grid[[currentBlockId(envir)]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gridApply() and blockApply()
###
### TODO: In theory, the best performance should be obtained when bplapply()
### uses a post office queue model. According to
### https://support.bioconductor.org/p/96856/#96888, this can be
### achieved by setting the nb of tasks to the nb of blocks (i.e. with
### BPPARAM=MulticoreParam(tasks=length(grid))). However, in practice, that
### seems to be slower than using tasks=0 (the default). Investigate this!
###

gridApply <- function(grid, FUN, ..., BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    if (!is(grid, "ArrayGrid"))
        stop(wmsg("'grid' must be an ArrayGrid object"))
    FUN <- match.fun(FUN)
    verbose <- normarg_verbose(verbose)

    FUN_WRAPPER <- function(bid, grid, verbose, FUN, ...)
    {
        if (verbose) {
            nblock <- length(grid)
            message("\\ Processing viewport ", bid, "/", nblock, " ... ",
                    appendLF=FALSE)
        }
        viewport <- grid[[bid]]
        set_grid_context(grid, bid, viewport)
        ans <- FUN(viewport, ...)
        if (verbose)
            message("OK")
        ans
    }
    S4Arrays:::bplapply2(seq_along(grid), FUN_WRAPPER, grid, verbose,
                                          FUN, ..., BPPARAM=BPPARAM)
}

viewportApply <- function(...)
{
    .Defunct("gridApply")
    gridApply(...)
}

blockApply <- function(x, FUN, ..., grid=NULL, as.sparse=FALSE,
                                    BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    FUN <- match.fun(FUN)
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be FALSE, TRUE, or NA"))
    verbose <- normarg_verbose(verbose)

    FUN_WRAPPER <- function(viewport,
                            FUN, x, as.sparse, verbose, verbose_read_block, ...)
    {
        effective_grid <- effectiveGrid()
        current_block_id <- currentBlockId()
        if (verbose) {
            x_is_sparse <- is_sparse(x)
            nblock <- length(effective_grid)
            block <- verbose_read_block(x, viewport, x_is_sparse,
                                        as.sparse, current_block_id, nblock)
        } else {
            block <- read_block(x, viewport, as.sparse=as.sparse)
        }
        set_grid_context(effective_grid, current_block_id, viewport)
        if (verbose)
            message("\\ Processing it ... ", appendLF=FALSE)
        ans <- FUN(block, ...)
        if (verbose)
            message("OK")
        ans
    }
    gridApply(grid, FUN_WRAPPER,
              FUN, x, as.sparse, verbose, verbose_read_block, ...,
              BPPARAM=BPPARAM, verbose=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gridReduce() and blockReduce()
###
### Two Reduce-like functions.
###

gridReduce <- function(FUN, grid, init, ..., BREAKIF=NULL, verbose=NA)
{
    FUN <- match.fun(FUN)
    if (!is(grid, "ArrayGrid"))
        stop(wmsg("'grid' must be an ArrayGrid object"))
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    verbose <- normarg_verbose(verbose)

    nblock <- length(grid)
    for (bid in seq_len(nblock)) {
        viewport <- grid[[bid]]
        set_grid_context(grid, bid, viewport)
        if (verbose)
            message("\\ Processing viewport ", bid, "/", nblock, " ... ",
                    appendLF=FALSE)
        init <- FUN(viewport, init, ...)
        if (verbose)
            message("OK")
        if (!is.null(BREAKIF) && BREAKIF(init)) {
            if (verbose)
                message("BREAK condition encountered")
            break
        }
    }
    init
}

viewportReduce <- function(...)
{
    .Defunct("gridReduce")
    gridReduce(...)
}

blockReduce <- function(FUN, x, init, ..., BREAKIF=NULL,
                        grid=NULL, as.sparse=FALSE, verbose=NA)
{
    FUN <- match.fun(FUN)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be FALSE, TRUE, or NA"))
    verbose <- normarg_verbose(verbose)

    FUN_WRAPPER <- function(viewport, init,
                            FUN, x, as.sparse, verbose, verbose_read_block, ...)
    {
        effective_grid <- effectiveGrid()
        current_block_id <- currentBlockId()
        if (verbose) {
            x_is_sparse <- is_sparse(x)
            nblock <- length(effective_grid)
            block <- verbose_read_block(x, viewport, x_is_sparse,
                                        as.sparse, current_block_id, nblock)
        } else {
            block <- read_block(x, viewport, as.sparse=as.sparse)
        }
        set_grid_context(effective_grid, current_block_id, viewport)
        if (verbose)
            message("\\ Processing it ... ", appendLF=FALSE)
        init <- FUN(block, init, ...)
        if (verbose)
            message("OK")
        init
    }
    if (!is.null(BREAKIF) && verbose) {
        BREAKIF_WRAPPER <- function(init)
        {
            ok <- BREAKIF(init)
            if (ok)
                message("BREAK condition encountered")
            ok
        }
    } else {
        BREAKIF_WRAPPER <- BREAKIF
        
    }
    gridReduce(FUN_WRAPPER, grid, init,
               FUN, x, as.sparse, verbose, verbose_read_block, ...,
               BREAKIF=BREAKIF_WRAPPER, verbose=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Two specialized internal helpers for the 2D case
###

### Return a list with one list element per row in the grid.
block_apply_by_grid_row <- function(INIT, INIT_args, FUN, FUN_args, x,
                                    grid=NULL, as.sparse=FALSE,
                                    BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    INIT <- match.fun(INIT)
    FUN <- match.fun(FUN)
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be a FALSE, TRUE, or NA"))
    verbose <- normarg_verbose(verbose)

    process_grid_row <- function(init, FUN, FUN_args, x,
                                 grid, i, as.sparse, verbose)
    {
        grid_nrow <- nrow(grid)
        grid_ncol <- ncol(grid)
        ## Inner loop on the grid columns. Sequential.
        for (j in seq_len(grid_ncol)) {
            viewport <- grid[[i, j]]
            block <- read_block(x, viewport, as.sparse=as.sparse)
            set_grid_context(grid, NULL, viewport)
            if (verbose)
                message("Processing block [[", i, "/", grid_nrow, ", ",
                                               j, "/", grid_ncol, "]] ... ",
                        appendLF=FALSE)
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

### Return a list with one list element per column in the grid.
block_apply_by_grid_col <- function(INIT, INIT_args, FUN, FUN_args, x,
                                    grid=NULL, as.sparse=FALSE,
                                    BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    INIT <- match.fun(INIT)
    FUN <- match.fun(FUN)
    grid <- normarg_grid(grid, x)
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be a FALSE, TRUE, or NA"))
    verbose <- normarg_verbose(verbose)

    process_grid_col <- function(init, FUN, FUN_args, x,
                                 grid, j, as.sparse, verbose)
    {
        grid_nrow <- nrow(grid)
        grid_ncol <- ncol(grid)
        ## Inner loop on the grid rows. Sequential.
        for (i in seq_len(grid_nrow)) {
            viewport <- grid[[i, j]]
            block <- read_block(x, viewport, as.sparse=as.sparse)
            set_grid_context(grid, NULL, viewport)
            if (verbose)
                message("Processing block [[", i, "/", grid_nrow, ", ",
                                               j, "/", grid_ncol, "]] ... ",
                        appendLF=FALSE)
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
### OLD - Walking on the blocks
### OLD -
### OLD - 3 utility functions to process array-like objects by block.
### OLD -
### OLD - Still used by the DelayedMatrixStats package.

### An lapply-like function.
block_APPLY <- function(x, APPLY, ..., sink=NULL, block_maxlen=NULL)
{
    APPLY <- match.fun(APPLY)
    x_dim <- dim(x)
    if (any(x_dim == 0L)) {
        chunk_grid <- NULL
    } else {
        ## Using chunks of length 1 (i.e. max resolution chunk grid) is just
        ## a trick to make sure that defaultAutoGrid() returns linear blocks.
        chunk_grid <- RegularArrayGrid(x_dim, rep.int(1L, length(x_dim)))
    }
    grid <- defaultAutoGrid(x, block_maxlen, chunk_grid,
                               block.shape="first-dim-grows-first")
    nblock <- length(grid)
    lapply(seq_len(nblock),
        function(bid) {
            if (get_verbose_block_processing())
                message("Processing block ", bid, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[bid]]
            block <- read_block(x, viewport)
            block_ans <- APPLY(block, ...)
            if (!is.null(sink)) {
                write_block(sink, viewport, block_ans)
                block_ans <- NULL
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        })
}

### A convenience wrapper around block_APPLY() to process a matrix-like
### object by block of columns.
colblock_APPLY <- function(x, APPLY, ..., sink=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    APPLY <- match.fun(APPLY)
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_maxlen <- max(getAutoBlockLength(type(x)), x_dim[[1L]])
    block_APPLY(x, APPLY, ..., sink=sink, block_maxlen=block_maxlen)
}

