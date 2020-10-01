### =========================================================================
### Block processing utilities
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
    set_user_option("auto.BPPARAM", BPPARAM)
}

getAutoBPPARAM <- function() get_user_option("auto.BPPARAM")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set/get grid context for the current block of a blockApply(),
### viewportApply(), or blockReduce() loop
###

### Exported (used in beachmat).
set_grid_context <- function(effective_grid, current_block_id,
                             envir=parent.frame(1))
{
    assign(".effective_grid", effective_grid, envir=envir)
    assign(".current_block_id", current_block_id, envir=envir)
}

.backward_compat <- function(envir, funname)
{
    if (!(is.array(envir) || is(envir, "SparseArraySeed")))
        return(envir)
    msg <- c("starting with DelayedArray 0.15.12, passing 'block' ",
             "to ", funname, "() is no longer needed and will soon ",
             "be considered an error")
    .Deprecated(msg=c("  ", wmsg(msg)))
    parent.frame(3)
}

.grid_context_not_found <- c(
    "Grid context not found for the current block. ",
    "Are we in a blockApply(), viewportApply(), or blockReduce() loop?"
)

.explain_proper_use <- function(funname)
    paste0("Note that ", funname, "() can only be called from **within** ",
           "the callback function of a blockApply() or viewportApply() ",
           "loop ('FUN' argument), or from the callback functions of a ",
           "blockReduce() loop ('FUN' and 'BREAKIF' arguments).")

.suggest_set_grid_context <- c(
    "If you need to be able to test/debug your callback function 'FUN' ",
    "(or 'BREAKIF') as a standalone function, call set_grid_context() to ",
    "set an arbitrary context **right before** calling the callback function."
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
### viewportApply() and blockApply()
###
### TODO: In theory, the best performance should be obtained when bplapply()
### uses a post office queue model. According to
### https://support.bioconductor.org/p/96856/#96888, this can be
### achieved by setting the nb of tasks to the nb of blocks (i.e. with
### BPPARAM=MulticoreParam(tasks=length(grid))). However, in practice, that
### seems to be slower than using tasks=0 (the default). Investigate this!
###

viewportApply <- function(grid, FUN, ..., BPPARAM=getAutoBPPARAM())
{
    if (!is(grid, "ArrayGrid"))
        stop(wmsg("'grid' must be an ArrayGrid object"))
    nblock <- length(grid)
    FUN <- match.fun(FUN)
    FUN0 <- function(bid, grid, nblock, FUN, ...) {
        if (get_verbose_block_processing()) {
            message("Processing block ", bid, "/", nblock, " ... ",
                    appendLF=FALSE)
            on.exit(message("OK"))
        }
        viewport <- grid[[bid]]
        set_grid_context(grid, bid)
        FUN(viewport, ...)
    }
    bplapply2(seq_len(nblock), FUN0, grid, nblock, FUN, ..., BPPARAM=BPPARAM)
}

blockApply <- function(x, FUN, ..., grid=NULL, as.sparse=FALSE,
                                    BPPARAM=getAutoBPPARAM())
{
    grid <- normarg_grid(grid, x)
    FUN <- match.fun(FUN)
    FUN0 <- function(viewport, x, as.sparse, FUN, ...) {
        block <- read_block(x, viewport, as.sparse=as.sparse)
        set_grid_context(effectiveGrid(), currentBlockId())
        FUN(block, ...)
    }
    viewportApply(grid, FUN0, x, as.sparse, FUN, ..., BPPARAM=BPPARAM)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### blockReduce()
###
### A Reduce-like function.
###

blockReduce <- function(FUN, x, init, BREAKIF=NULL, grid=NULL, as.sparse=FALSE)
{
    FUN <- match.fun(FUN)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    grid <- normarg_grid(grid, x)
    nblock <- length(grid)
    for (bid in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", bid, "/", nblock, " ... ",
                    appendLF=FALSE)
        viewport <- grid[[bid]]
        block <- read_block(x, viewport, as.sparse=as.sparse)
        set_grid_context(grid, bid)
        init <- FUN(block, init)
        if (get_verbose_block_processing())
            message("OK")
        if (!is.null(BREAKIF) && BREAKIF(init)) {
            if (get_verbose_block_processing())
                message("BREAK condition encountered")
            break
        }
    }
    init
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

