### =========================================================================
### Block processing utilities
### -------------------------------------------------------------------------
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

.normarg_chunk.grid <- function(chunk.grid, x)
{
    if (is.null(chunk.grid))
        return(chunkGrid(x))
    if (!is(chunk.grid, "ArrayGrid"))
        stop(wmsg("'chunk.grid' must be NULL or an ArrayGrid object"))
    if (!identical(refdim(chunk.grid), dim(x)))
        stop(wmsg("'chunk.grid' is incompatible with 'x'"))
    chunk.grid
}

### Return an "optimal" grid for block processing.
### The grid is returned as an ArrayGrid object on reference array 'x'.
### The grid elements define the blocks that will be used for processing 'x'
### by block. The grid is "optimal" in the sense that:
###  - It's "compatible" with the chunk grid (i.e. with 'chunkGrid(x)' or
###    with the chunk grid supplied via the 'chunk.grid' argument), that is,
###    the chunks are contained in the blocks. Said otherwise, chunks never
###    cross block boundaries.
###  - Its "resolution" is such that the blocks have a length that is as
###    close as possibe to (but does not exceed) 'block.maxlength'.
###    An exception is when some chunks are already >= 'block.maxlength',
###    in which case the returned grid is the same as the chunk grid.
### Note that the returned grid is regular (i.e. RegularArrayGrid object)
### unless the chunk grid is not regular (i.e. is an ArbitraryArrayGrid
### object).
defaultGrid <- function(x, block.maxlength=NULL,
                           chunk.grid=NULL,
                           block.shape=c("hypercube", "linear"))
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must have dimensions"))
    block_maxlen <- .normarg_block.maxlength(block.maxlength, type(x))
    chunk_grid <- .normarg_chunk.grid(chunk.grid, x)
    block_shape <- match.arg(block.shape)
    ## If 'x' is empty, we return a grid with a single (empty) block that
    ## has the dimensions of 'x'.
    if (any(x_dim == 0L))
        return(RegularArrayGrid(x_dim))
    if (is.null(chunk_grid)) {
        ans <- make_RegularArrayGrid_of_capped_length_blocks(
                           x_dim, block_maxlen, block_shape=block_shape)
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

.as_array_or_matrix <- function(x)
{
    if (length(dim(x)) == 2L)
        return(as.matrix(x))
    as.array(x)
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
### block_which()
###

### Return a numeric matrix like one returned by base::arrayInd(), that is,
### a matrix where each row is an n-uplet representing an array index.
### 'x' is **trusted** to be a logical array-like object.
block_which <- function(x, grid=NULL)
{
    FUN <- function(block) {
        b <- currentBlockId(block)
        m <- base::which(block)
        mapToRef(rep.int(b, length(m)), m, effectiveGrid(block), linear=TRUE)
    }
    block_results <- blockApply(x, FUN, grid=grid)
    aind <- do.call(rbind, block_results)
    aind_as_list <- lapply(ncol(aind):1, function(j) aind[ , j])
    oo <- do.call(order, aind_as_list)
    aind[oo, , drop=FALSE]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### block_extract_array_elements()
###

### Linear single bracket subsetting (e.g. x[5:2]) of an array-like object.
### Return an atomic vector.
### 'x' is **trusted** to be an array-like object.
### 'i' is **trusted** to be an integer vector representing a linear index
### of valid positions in 'x'.
block_extract_array_elements <- function(x, i, grid=NULL)
{
    i_len <- length(i)
    if (i_len == 0L)
        return(as.vector(extract_empty_array(x)))
    if (i_len == 1L)
        return(extract_array_element(x, i))

    ## We don't want to use blockApply() here because it would walk on the
    ## entire grid of blocks, which is not necessary. We only need to walk
    ## on the blocks touched by linear index 'i', that is, on the blocks
    ## that contain array elements located at the positions corresponding
    ## to linear index 'i'.
    grid <- .normarg_grid(grid, x)
    nblock <- length(grid)
    x_dim <- dim(x)
    majmin <- mapToGrid(arrayInd(i, x_dim), grid, linear=TRUE)
    minor_by_block <- split(majmin$minor, majmin$major)
    res <- lapply(seq_along(minor_by_block),
        function(k) {
            b <- as.integer(names(minor_by_block)[[k]])
            m <- minor_by_block[[k]]
            if (get_verbose_block_processing())
                message("Visiting block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            ## We don't need to load the entire block if there is only 1
            ## value to extract from it.
            if (length(m) == 1L) {
                i2 <- linearInd(mapToRef(b, m, grid, linear=TRUE), x_dim)
                block_ans <- extract_array_element(x, i2)
            } else {
                block <- extract_block(x, grid[[b]])
                if (!is.array(block))
                    block <- as.array(block)
                block_ans <- block[m]
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
    })
    unsplit(res, majmin$major)
}


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
        ## a trick to make sure that defaultGrid() returns linear blocks.
        chunk_grid <- RegularArrayGrid(x_dim, rep.int(1L, length(x_dim)))
    }
    grid <- defaultGrid(x, block_maxlen, chunk_grid, block.shape="linear")
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

### Split the array-like object into blocks, then realize and write one block
### at a time to disk.
write_array_to_sink <- function(x, sink)
{
    stopifnot(identical(dim(x), dim(sink)))
    block_APPLY(DelayedArray(x), identity, sink=sink)
}

