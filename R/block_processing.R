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
get_block_length <- function(type)
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
block_APPLY <- function(x, APPLY, ..., sink=NULL, block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    if (is.null(block_len))
        block_len <- get_block_length(type(x))
    blocks <- ArrayBlocks(dim(x), block_len)
    nblock <- length(blocks)
    lapply(seq_len(nblock),
        function(i) {
            if (get_verbose_block_processing())
                message("Processing block ", i, "/", nblock, " ... ",
                        appendLF=FALSE)
            block_ranges <- get_block_ranges(blocks, i)
            Nindex <- make_Nindex_from_block_ranges(block_ranges, blocks@dim)
            subarray <- subset_by_Nindex(x, Nindex)
            if (!is.array(subarray))
                subarray <- .as_array_or_matrix(subarray)
            block_ans <- APPLY(subarray, ...)
            if (!is.null(sink)) {
                write_to_sink(block_ans, sink, offsets=start(block_ranges))
                block_ans <- NULL
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        })
}

### A mapply-like function for conformable arrays.
block_MAPPLY <- function(MAPPLY, ..., sink=NULL, block_len=NULL)
{
    MAPPLY <- match.fun(MAPPLY)
    dots <- unname(list(...))
    dims <- sapply(dots, dim)
    x_dim <- dims[ , 1L]
    if (!all(dims == x_dim))
        stop("non-conformable arrays")
    if (is.null(block_len)) {
        types <- unlist(lapply(dots, type))
        block_len <- min(get_block_length(types))
    }
    blocks <- ArrayBlocks(x_dim, block_len)
    nblock <- length(blocks)
    lapply(seq_len(nblock),
        function(i) {
            if (get_verbose_block_processing())
                message("Processing block ", i, "/", nblock, " ... ",
                        appendLF=FALSE)
            block_ranges <- get_block_ranges(blocks, i)
            Nindex <- make_Nindex_from_block_ranges(block_ranges, blocks@dim)
            subarrays <- lapply(dots,
                function(x) {
                    subarray <- subset_by_Nindex(x, Nindex)
                    if (!is.array(subarray))
                        subarray <- .as_array_or_matrix(subarray)
                    subarray
                })
            block_ans <- do.call(MAPPLY, subarrays)
            if (!is.null(sink)) {
                write_to_sink(block_ans, sink, offsets=start(block_ranges))
                block_ans <- NULL
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
        })
}

### A Reduce-like function.
block_APPLY_and_COMBINE <- function(x, APPLY, COMBINE, init,
                                       BREAKIF=NULL, block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    COMBINE <- match.fun(COMBINE)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    if (is.null(block_len))
        block_len <- get_block_length(type(x))
    blocks <- ArrayBlocks(dim(x), block_len)
    nblock <- length(blocks)
    for (i in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", i, "/", nblock, " ... ",
                    appendLF=FALSE)
        subarray <- extract_array_block(x, blocks, i)
        if (!is.array(subarray))
            subarray <- .as_array_or_matrix(subarray)
        reduced <- APPLY(subarray)
        init <- COMBINE(i, subarray, init, reduced)
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
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    block_APPLY(x, APPLY, ..., sink=sink, block_len=block_len)
}

colblock_APPLY_and_COMBINE <- function(x, APPLY, COMBINE, init)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    block_APPLY_and_COMBINE(x, APPLY, COMBINE, init, block_len=block_len)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Block by block realization of an array-like object
###

### Exported!
setMethod("write_to_sink", c("DelayedArray", "RealizationSink"),
    function(x, sink, offsets=NULL)
    {
        if (!is.null(offsets))
            stop(wmsg("'offsets' must be NULL when the object to write ",
                      "is a DelayedArray object"))
        ## Semantically equivalent to 'write_to_sink(as.array(x), sink)'
        ## but uses block-processing so the full DelayedArray object is
        ## not realized at once in memory. Instead the object is first
        ## split into blocks and the blocks are realized and written to
        ## disk one at a time.
        block_APPLY(x, identity, sink=sink)
    }
)

### Exported!
setMethod("write_to_sink", c("ANY", "RealizationSink"),
    function(x, sink, offsets=NULL)
    {
        x <- as(x, "DelayedArray")
        write_to_sink(x, sink, offsets=offsets)
    }
)

