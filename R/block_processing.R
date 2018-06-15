### =========================================================================
### Internal utilities for block processing an array
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


### Default block size in bytes.
DEFAULT_BLOCK_SIZE <- 45000000L  # 45 Mb

### Atomic type sizes in bytes.
.TYPE_SIZES <- c(
    logical=4L,
    integer=4L,
    numeric=8L,
    double=8L,
    complex=16L,
    character=8L,  # just the overhead of a CHARSXP; doesn't account for the
                   # string data itself
    raw=1L
)

### NOT exported but used in HDF5Array!
get_default_block_maxlength <- function(type)
{
    type_size <- .TYPE_SIZES[type]
    idx <- which(is.na(type_size))
    if (length(idx) != 0L) {
        unsupported_types <- unique(type[idx])
        in1string <- paste0(unsupported_types, collapse=", ")
        stop("unsupported type(s): ",  in1string)
    }
    block_size <- getOption("DelayedArray.block.size",
                            default=DEFAULT_BLOCK_SIZE)
    if (!isSingleNumber(block_size) || block_size < 1)
        stop(wmsg("global option DelayedArray.block.size must be a ",
                  "single number >= 1"))
    ans <- block_size / type_size
    if (ans > .Machine$integer.max)
        stop(wmsg("Default block length is too big. Blocks of ",
                  "length > .Machine$integer.max are not supported yet. ",
                  "Please reduce the default block length by setting global ",
                  "option DelayedArray.block.size to a smaller value."))
    max(as.integer(ans), 1L)
}

### Guaranteed to return an integer >= 1.
.normarg_block.maxlength <- function(block.maxlength, type)
{
    if (is.null(block.maxlength))
        return(get_default_block_maxlength(type))
    if (!isSingleNumber(block.maxlength))
        stop(wmsg("'block.maxlength' must be a single integer or NULL"))
    if (block.maxlength < 1)
        stop(wmsg("'block.maxlength' cannot be < 1"))
    if (block.maxlength > .Machine$integer.max)
        stop(wmsg("'block.maxlength' is too big. Blocks of ",
                  "length > .Machine$integer.max are not supported yet. ",
                  "Please specify a smaller 'block.maxlength'."))
    as.integer(block.maxlength)
}

### Return a regular grid with linear blocks of length as close as possibe to
### (but not bigger than) 'block.maxlength'.
### TODO: Add 'block.shape' argument.
defaultGrid <- function(x, block.maxlength=NULL)
#                           block.shape=c("hypercube", "linear"))
{
    block_maxlen <- .normarg_block.maxlength(block.maxlength, type(x))
    block_shape <- "linear"
    chunk_grid <- NULL
    #block_shape <- match.arg(block.shape)
    #chunk_grid <- chunkGrid(x)
    if (is.null(chunk_grid)) {
        ans <- make_RegularArrayGrid_of_capped_length_blocks(
                           dim(x), block_maxlen, block_shape=block_shape)
        return(ans)
    }
    chunks_per_block <- max(block_maxlen %/% maxlength(chunk_grid), 1L)
    ratio <- get_spacings_for_capped_length_blocks(
                 dim(chunk_grid), chunks_per_block, block_shape=block_shape)
    downsample(chunk_grid, ratio)
}

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
### Walking on the blocks
###
### 2 utility functions to process array-like objects by block.
###

.as_array_or_matrix <- function(x)
{
    if (length(dim(x)) == 2L)
        return(as.matrix(x))
    as.array(x)
}

.normarg_grid <- function(grid, x)
{
    if (is.null(grid))
        return(defaultGrid(x))
    if (!is(grid, "ArrayGrid"))
        stop(wmsg("'grid' must be NULL or an ArrayGrid object"))
    if (!identical(refdim(grid), dim(x)))
        stop(wmsg("'grid' is incompatible with 'x'"))
    grid
}

### 'x' must be an array-like object.
### 'FUN' is the callback function to be applied to each block of array-like
### object 'x'. It must take at least 1 argument which is the current array
### block as an ordinary array or matrix.
### 'grid' must be an ArrayGrid object describing the block partitioning
### of 'x'. If not supplied, the grid returned by 'defaultGrid(x)' is used.
### The effective grid (i.e. 'grid' or 'defaultGrid(x)'), current block number,
### and current viewport (i.e. the ArrayViewport object describing the position
### of the current block w.r.t. the effective grid), can be obtained from
### within 'FUN' with 'effectiveGrid(block)', 'currentBlockId(block)', and
### 'currentViewport(block)', respectively.
### 'BPREDO' and 'BPPARAM' are passed to bplapply(). In theory, the best
### performance should be obtained when bplapply() uses a post office queue
### model. According to https://support.bioconductor.org/p/96856/#96888, this
### can be achieved by setting the nb of tasks to the nb of blocks (i.e. with
### BPPARAM=MulticoreParam(tasks=length(grid))). However, in practice, that
### seems to be slower than using tasks=0 (the default). Investigate this!
blockApply <- function(x, FUN, ..., grid=NULL, BPREDO=list(), BPPARAM=bpparam())
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
            block <- extract_block(x, viewport)
            if (!is.array(block))
                block <- .as_array_or_matrix(block)
            attr(block, "from_grid") <- grid
            attr(block, "block_id") <- b
            block_ans <- FUN(block, ...)
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        },
        BPREDO=BPREDO,
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
        block <- extract_block(x, viewport)
        if (!is.array(block))
            block <- .as_array_or_matrix(block)
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

### An lapply-like function.
block_APPLY <- function(x, APPLY, ..., sink=NULL, block_maxlen=NULL)
{
    APPLY <- match.fun(APPLY)
    grid <- defaultGrid(x, block_maxlen)
    nblock <- length(grid)
    lapply(seq_len(nblock),
        function(b) {
            if (get_verbose_block_processing())
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[b]]
            block <- extract_block(x, viewport)
            if (!is.array(block))
                block <- .as_array_or_matrix(block)
            block_ans <- APPLY(block, ...)
            if (!is.null(sink)) {
                write_block_to_sink(block_ans, sink, viewport)
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
    block_maxlen <- max(get_default_block_maxlength(type(x)), x_dim[[1L]])
    block_APPLY(x, APPLY, ..., sink=sink, block_maxlen=block_maxlen)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Block by block realization of an array-like object
###

### Exported!
### Split the array-like object into blocks, then realize and write one block
### at a time to disk.
write_array_to_sink <- function(x, sink)
{
    stopifnot(identical(dim(x), dim(sink)))
    block_APPLY(DelayedArray(x), identity, sink=sink)
}

