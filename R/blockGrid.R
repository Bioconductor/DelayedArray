### =========================================================================
### blockGrid() and family
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
###    the chunks are contained in the blocks. In other words, chunks never
###    cross block boundaries.
###  - Its "resolution" is such that the blocks have a length that is as
###    close as possibe to (but does not exceed) 'block.maxlength'.
###    An exception is when some chunks are already >= 'block.maxlength',
###    in which case the returned grid is the same as the chunk grid.
### Note that the returned grid is regular (i.e. RegularArrayGrid object)
### unless the chunk grid is not regular (i.e. is an ArbitraryArrayGrid
### object).
blockGrid <- function(x, block.maxlength=NULL,
                         chunk.grid=NULL,
                         block.shape=c("hypercube",
                                       "scale",
                                       "first-dim-grows-first",
                                       "last-dim-grows-first"))
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
        ans <- makeRegularArrayGridOfCappedLengthViewports(x_dim,
                                                           block_maxlen,
                                                           block_shape)
        return(ans)
    }
    chunks_per_block <- max(block_maxlen %/% maxlength(chunk_grid), 1L)
    ratio <- makeCappedVolumeBox(chunks_per_block, dim(chunk_grid), block_shape)
    downsample(chunk_grid, ratio)
}

