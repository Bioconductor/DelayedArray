### =========================================================================
### Internal utilities for block processing an array
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


### Default block size in bytes.
DEFAULT_BLOCK_SIZE <- 4500000L  # 4.5 Mb

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

### Used in HDF5Array!
get_max_block_length <- function(type)
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
    as.integer(block_size / type_size)
}

### Return a regular grid with linear blocks of length as close as possibe to
### (but not bigger than) 'max.block.length'.
### TODO: Add 'block.shape' argument.
defaultGrid <- function(x, max.block.length=NULL)
{
    if (is.null(max.block.length)) {
        max_block_len <- get_max_block_length(type(x))
    } else {
        if (!isSingleNumber(max.block.length))
            stop("'max.block.length' must be a single integer or NULL")
        max_block_len <- as.integer(max.block.length)
    }
    make_RegularArrayGrid_of_capped_length_blocks(dim(x), max_block_len,
                                                  block_shape="linear")
}

### Used in HDF5Array!
get_verbose_block_processing <- function()
{
    getOption("DelayedArray.verbose.block.processing", default=FALSE)
}

### Used in HDF5Array!
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
### 3 utility functions to process array-like objects by block.
###

.as_array_or_matrix <- function(x)
{
    if (length(dim(x)) == 2L)
        return(as.matrix(x))
    as.array(x)
}

### An lapply-like function.
block_APPLY <- function(x, APPLY, ..., sink=NULL, max_block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    grid <- defaultGrid(x, max_block_len)
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
### block as an ordinay array or matrix.
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

### A mapply-like function for conformable arrays.
block_MAPPLY <- function(MAPPLY, ..., sink=NULL, max_block_len=NULL)
{
    MAPPLY <- match.fun(MAPPLY)
    dots <- unname(list(...))
    dims <- vapply(dots, dim, integer(2))
    if (!all(dims == dims[ , 1L]))
        stop("non-conformable arrays")
    if (is.null(max_block_len)) {
        types <- unlist(lapply(dots, type))
        max_block_len <- min(get_max_block_length(types))
    }
    grid <- defaultGrid(x, max_block_len)
    nblock <- length(grid)
    lapply(seq_len(nblock),
        function(b) {
            if (get_verbose_block_processing())
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[b]]
            blocks <- lapply(dots,
                function(x) {
                    block <- extract_block(x, viewport)
                    if (!is.array(block))
                        block <- .as_array_or_matrix(block)
                    block
                })
            block_ans <- do.call(MAPPLY, blocks)
            if (!is.null(sink)) {
                write_block_to_sink(block_ans, sink, viewport)
                block_ans <- NULL
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        })
}

### A Reduce-like function.
block_APPLY_and_COMBINE <- function(x, APPLY, COMBINE, init,
                                       BREAKIF=NULL, max_block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    COMBINE <- match.fun(COMBINE)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    grid <- defaultGrid(x, max_block_len)
    nblock <- length(grid)
    for (b in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF=FALSE)
        block <- extract_block(x, grid[[b]])
        if (!is.array(block))
            block <- .as_array_or_matrix(block)
        reduced <- APPLY(block)
        init <- COMBINE(b, block, init, reduced)
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
### Walking on the blocks of columns
###
### 2 convenience wrappers around block_APPLY() and block_APPLY_and_COMBINE()
### to process a matrix-like object by block of columns.
###

colblock_APPLY <- function(x, APPLY, ..., sink=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    APPLY <- match.fun(APPLY)
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    max_block_len <- max(get_max_block_length(type(x)), x_dim[[1L]])
    block_APPLY(x, APPLY, ..., sink=sink, max_block_len=max_block_len)
}

colblock_APPLY_and_COMBINE <- function(x, APPLY, COMBINE, init)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    max_block_len <- max(get_max_block_length(type(x)), x_dim[[1L]])
    block_APPLY_and_COMBINE(x, APPLY, COMBINE, init,
                            max_block_len=max_block_len)
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

