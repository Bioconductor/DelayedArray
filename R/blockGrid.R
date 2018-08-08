### =========================================================================
### blockGrid() and family
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get/setDefaultBlockSize()
###
### The default block size must be specified in bytes.
###

getDefaultBlockSize <- function()
{
    size <- getOption("DelayedArray.block.size")
    if (!isSingleNumber(size) || size < 1)
        stop(wmsg("Global option DelayedArray.block.size ",
                  "should be a single number >= 1. ",
                  "Fix it with setDefaultBlockSize()."))
    size
}

### We set the default block size to 100Mb by default.
setDefaultBlockSize <- function(size=1e8)
{
    if (!isSingleNumber(size) || size < 1)
        stop(wmsg("the block size must be a single number >= 1"))
    prev_size <- getOption("DelayedArray.block.size")
    if (isSingleNumber(prev_size)) {
        previous_was <- c(" (was ", prev_size, ")")
    } else {
        previous_was <- ""
    }
    options(DelayedArray.block.size=size)
    message("default block size set to ", size, " bytes", previous_was)
    invisible(size)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getDefaultBlockLength()
###

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

getDefaultBlockLength <- function(type)
{
    type_size <- .TYPE_SIZES[type]
    idx <- which(is.na(type_size))
    if (length(idx) != 0L) {
        unsupported_types <- unique(type[idx])
        in1string <- paste0(unsupported_types, collapse=", ")
        stop("unsupported type(s): ",  in1string)
    }
    block_size <- getDefaultBlockSize()
    ans <- block_size / type_size
    if (ans > .Machine$integer.max)
        stop(wmsg("Default block length is too big. Blocks of ",
                  "length > .Machine$integer.max are not supported yet. ",
                  "Please reduce the default block length by reducing the ",
                  "default block size with setDefaultBlockSize()."))
    max(as.integer(ans), 1L)
}

### Guaranteed to return an integer >= 1.
.normarg_block.length <- function(block.length, type)
{
    if (is.null(block.length))
        return(getDefaultBlockLength(type))
    if (!isSingleNumber(block.length))
        stop(wmsg("'block.length' must be a single integer or NULL"))
    if (block.length < 1)
        stop(wmsg("'block.length' must be >= 1"))
    if (block.length > .Machine$integer.max)
        stop(wmsg("'block.length' is too big. Blocks of ",
                  "length > .Machine$integer.max are not supported yet. ",
                  "Please specify a smaller 'block.length'."))
    as.integer(block.length)
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

### Return a grid that is "optimal" for block processing of array-like
### object 'x'.
### The grid is returned as an ArrayGrid object on reference array 'x'.
### The grid elements define the blocks that will be used for processing 'x'
### by block. The grid is "optimal" in the sense that:
###  - It's "compatible" with the chunk grid (i.e. with 'chunkGrid(x)' or
###    with the chunk grid supplied via the 'chunk.grid' argument), that is,
###    the chunks are contained in the blocks. In other words, chunks never
###    cross block boundaries.
###  - Its "resolution" is such that the blocks have a length that is as
###    close as possibe to (but does not exceed) 'block.length'.
###    An exception is when some chunks are already >= 'block.length',
###    in which case the returned grid is the same as the chunk grid.
### Note that the returned grid is regular (i.e. RegularArrayGrid object)
### unless the chunk grid is not regular (i.e. is an ArbitraryArrayGrid
### object).
blockGrid <- function(x, block.length=NULL,
                         chunk.grid=NULL,
                         block.shape=c("hypercube",
                                       "scale",
                                       "first-dim-grows-first",
                                       "last-dim-grows-first"))
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must have dimensions"))
    block_len <- .normarg_block.length(block.length, type(x))
    chunk_grid <- .normarg_chunk.grid(chunk.grid, x)
    block_shape <- match.arg(block.shape)
    ## If 'x' is empty, we return a grid with a single (empty) block that
    ## has the dimensions of 'x'.
    if (any(x_dim == 0L))
        return(RegularArrayGrid(x_dim))
    if (is.null(chunk_grid)) {
        ans <- makeRegularArrayGridOfCappedLengthViewports(x_dim,
                                                           block_len,
                                                           block_shape)
        return(ans)
    }
    chunks_per_block <- max(block_len %/% maxlength(chunk_grid), 1L)
    ratio <- makeCappedVolumeBox(chunks_per_block, dim(chunk_grid), block_shape)
    downsample(chunk_grid, ratio)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Two additional functions specific to the 2-dimensional case
###
### Both return a RegularArrayGrid object.
###

.get_default_nrow <- function(x_dim, block.length, x_type)
{
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    block_len <- .normarg_block.length(block.length, x_type)
    nrow <- block_len %/% x_ncol
    if (nrow < 1L)
        return(1L)
    if (nrow > x_nrow)
        return(x_nrow)
    nrow
}

### Define blocks of full rows.
rowGrid <- function(x, nrow=NULL, block.length=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (is.null(nrow)) {
        nrow <- .get_default_nrow(x_dim, block.length, type(x))
        spacings <- c(nrow, x_ncol)
    } else {
        if (!is.null(block.length))
            warning("'block.length' is ignored when 'nrow' is not NULL")
        if (!isSingleNumber(nrow))
            stop(wmsg("'nrow' must be a single integer or NULL"))
        nrow <- as.integer(nrow)
        if (nrow < 1L || nrow > x_nrow)
            stop(wmsg("'nrow' must be >= 1 and <= nrow(x)"))
        spacings <- c(nrow, x_ncol)
        if (prod(spacings) > .Machine$integer.max)
            stop(wmsg("'nrow' is too big. Blocks of length > ",
                      ".Machine$integer.max are not supported yet. ",
                      "Please specify a smaller 'nrow'."))
    }
    RegularArrayGrid(x_dim, spacings)
}

### Define blocks of full columns.
colGrid <- function(x, ncol=NULL, block.length=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (is.null(ncol)) {
        ncol <- .get_default_nrow(rev(x_dim), block.length, type(x))
        spacings <- c(x_nrow, ncol)
    } else {
        if (!is.null(block.length))
            warning("'block.length' is ignored when 'ncol' is not NULL")
        if (!isSingleNumber(ncol))
            stop(wmsg("'ncol' must be a single integer or NULL"))
        ncol <- as.integer(ncol)
        if (ncol < 1L || ncol > x_ncol)
            stop(wmsg("'ncol' must be >= 1 and <= ncol(x)"))
        spacings <- c(x_nrow, ncol)
        if (prod(spacings) > .Machine$integer.max)
            stop(wmsg("'ncol' is too big. Blocks of length > ",
                      ".Machine$integer.max are not supported yet. ",
                      "Please specify a smaller 'ncol'."))
    }
    RegularArrayGrid(x_dim, spacings)
}

