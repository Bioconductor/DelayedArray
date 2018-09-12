### =========================================================================
### Block processing utilities
### -------------------------------------------------------------------------
###


### NOT exported but used in HDF5Array!
get_verbose_block_processing <- function()
{
    getOption("DelayedArray.verbose.block.processing", default=FALSE)
}

### NOT exported but used in HDF5Array!
set_verbose_block_processing <- function(verbose)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    old_verbose <- get_verbose_block_processing()
    options(DelayedArray.verbose.block.processing=verbose)
    old_verbose
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoGridMaker()
###

### We set the automatic grid maker to blockGrid() by default.
setAutoGridMaker <- function(GRIDMAKER="blockGrid")
{
    match.fun(GRIDMAKER)  # sanity check
    set_user_option("auto.grid.maker", GRIDMAKER)
}

getAutoGridMaker <- function() get_user_option("auto.grid.maker")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoBPPARAM()
###

### By default (i.e. when no argument is specified), we set the automatic
### BPPARAM to SerialParam() on Windows and to MulticoreParam() on other
### platforms.
### The reason we use SerialParam() instead of SnowParam() on Windows is
### that the latter introduces **a lot** of overhead in the context of block
### processing. See https://github.com/Bioconductor/BiocParallel/issues/78
### To the point that disabling parallel evaluation (by using SerialParam())
### is still much faster than parallel evaluation with SnowParam().
setAutoBPPARAM <- function(BPPARAM=NULL)
{
    if (is.null(BPPARAM)) {
        if (.Platform$OS.type == "windows") {
            BPPARAM <- BiocParallel::SerialParam()
        } else {
            BPPARAM <- BiocParallel::MulticoreParam()
        }
    } else {
        if (!is(BPPARAM, "BiocParallelParam"))
            stop(wmsg("'BPPARAM' must be a BiocParallelParam object"))
    }
    set_user_option("auto.BPPARAM", BPPARAM)
}

getAutoBPPARAM <- function() get_user_option("auto.BPPARAM")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks
###
### 2 utility functions to process array-like objects by block.
###

.normarg_grid <- function(grid, x)
{
    if (is.null(grid)) {
        etc <- c("Please use setAutoGridMaker() ",
                 "to set a valid automatic grid maker.")
        GRIDMAKER <- match.fun(getAutoGridMaker())
        grid <- try(GRIDMAKER(x), silent=TRUE)
        if (is(grid, "try-error"))
            stop(wmsg("The current automatic grid maker returned an ",
                      "error when called on 'x'. ", etc))
        if (!is(grid, "ArrayGrid"))
            stop(wmsg("The current automatic grid maker didn't return an ",
                      "ArrayGrid object. ", etc))
        if (!identical(refdim(grid), dim(x)))
            stop(wmsg("The current automatic grid maker returned a grid ",
                      "that is incompatible with 'x'. ", etc))
    } else {
        if (!is(grid, "ArrayGrid"))
            stop(wmsg("'grid' must be NULL or an ArrayGrid object"))
        if (!identical(refdim(grid), dim(x)))
            stop(wmsg("'grid' is incompatible with 'x'"))
    }
    grid
}

### 'x' must be an array-like object.
### 'FUN' is the callback function to be applied to each block of array-like
### object 'x'. It must take at least 1 argument which is the current array
### block as an ordinary array or matrix.
### 'grid' must be an ArrayGrid object describing the block partitioning
### of 'x'. If not supplied, the grid returned by 'blockGrid(x)' is used.
### The effective grid (i.e. 'grid' or 'blockGrid(x)'), current block number,
### and current viewport (i.e. the ArrayViewport object describing the position
### of the current block w.r.t. the effective grid), can be obtained from
### within 'FUN' with 'effectiveGrid(block)', 'currentBlockId(block)', and
### 'currentViewport(block)', respectively.
### 'BPPARAM' is passed to bplapply(). In theory, the best performance should
### be obtained when bplapply() uses a post office queue model. According to
### https://support.bioconductor.org/p/96856/#96888, this can be achieved by
### setting the nb of tasks to the nb of blocks (i.e. with
### BPPARAM=MulticoreParam(tasks=length(grid))). However, in practice, that
### seems to be slower than using tasks=0 (the default). Investigate this!
blockApply <- function(x, FUN, ..., grid=NULL, BPPARAM=getAutoBPPARAM())
{
    FUN <- match.fun(FUN)
    grid <- .normarg_grid(grid, x)
    nblock <- length(grid)
    bplapply(seq_len(nblock),
        function(b) {
            if (get_verbose_block_processing())
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[b]]
            block <- read_block(x, viewport)
            attr(block, "from_grid") <- grid
            attr(block, "block_id") <- b
            block_ans <- FUN(block, ...)
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        },
        BPPARAM=BPPARAM
    )
}

### A Reduce-like function. Not parallelized yet.
blockReduce <- function(FUN, x, init, BREAKIF=NULL, grid=NULL)
{
    FUN <- match.fun(FUN)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    grid <- .normarg_grid(grid, x)
    nblock <- length(grid)
    for (b in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF=FALSE)
        viewport <- grid[[b]]
        block <- read_block(x, viewport)
        attr(block, "from_grid") <- grid
        attr(block, "block_id") <- b
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

effectiveGrid <- function(block)
{
    if (!is.array(block))
        stop("'block' must be an ordinary array")
    if (!("from_grid" %in% names(attributes(block))))
        stop(wmsg("'block' has no \"from_grid\" attribute. ",
                  "Was effectiveGrid() called in a blockApply() loop?"))
    attr(block, "from_grid", exact=TRUE)
}

currentBlockId <- function(block)
{
    if (!is.array(block))
        stop("'block' must be an ordinary array")
    if (!("block_id" %in% names(attributes(block))))
        stop(wmsg("'block' has no \"block_id\" attribute. ",
                  "Was currentBlockId() called in a blockApply() loop?"))
    attr(block, "block_id", exact=TRUE)
}

currentViewport <- function(block)
    effectiveGrid(block)[[currentBlockId(block)]]


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OLD - Walking on the blocks
### OLD -
### OLD - 3 utility functions to process array-like objects by block.
### OLD -
### OLD - Still used by write_array_to_sink() below and by the
### OLD - DelayedMatrixStats package.

### An lapply-like function.
block_APPLY <- function(x, APPLY, ..., sink=NULL, block_maxlen=NULL)
{
    APPLY <- match.fun(APPLY)
    x_dim <- dim(x)
    if (any(x_dim == 0L)) {
        chunk_grid <- NULL
    } else {
        ## Using chunks of length 1 (i.e. max resolution chunk grid) is just
        ## a trick to make sure that blockGrid() returns linear blocks.
        chunk_grid <- RegularArrayGrid(x_dim, rep.int(1L, length(x_dim)))
    }
    grid <- blockGrid(x, block_maxlen, chunk_grid,
                         block.shape="first-dim-grows-first")
    nblock <- length(grid)
    lapply(seq_len(nblock),
        function(b) {
            if (get_verbose_block_processing())
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[b]]
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

