### =========================================================================
### Internal utilities for block processing an array
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayBlocks objects
###

setClass("ArrayBlocks",
    representation(
        dim="integer",
        N="integer",
        by="integer"
    )
)

### Return an ArrayBlocks object i.e. a collection of subarrays of the
### original array with the following properties:
###   (a) The collection of blocks is a partitioning of the original array
###       i.e. the blocks fully cover it and don't overlap each other.
###   (b) Each block is made of adjacent elements in the original array.
###   (c) Each block has a length (i.e. nb of elements) <= 'max_block_len'.
ArrayBlocks <- function(dim, max_block_len)
{
    p <- cumprod(dim)
    w <- which(p <= max_block_len)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N > length(dim)) {
        by <- 1L
    } else if (N == 1L) {
        by <- max_block_len
    } else {
        by <- max_block_len %/% as.integer(p[[N - 1L]])
    }
    new("ArrayBlocks", dim=dim, N=N, by=by)
}

.get_ArrayBlocks_inner_length <- function(x)
{
    ndim <- length(x@dim)
    if (x@N > ndim)
        return(if (any(x@dim == 0L)) 0L else 1L)
    inner_len <- x@dim[[x@N]] %/% x@by
    last_inner_block_len <- x@dim[[x@N]] %% x@by
    if (last_inner_block_len != 0L)
        inner_len <- inner_len + 1L
    inner_len
}

.get_ArrayBlocks_outer_length <- function(x)
{
    ndim <- length(x@dim)
    if (x@N >= ndim)
        return(1L)
    outer_dim <- x@dim[(x@N + 1L):ndim]
    prod(outer_dim)
}

### Return the number of blocks in 'x'.
setMethod("length", "ArrayBlocks",
    function(x)
        .get_ArrayBlocks_inner_length(x) * .get_ArrayBlocks_outer_length(x)
)

### NOT exported but used in HDF5Array!
### Return a "multidimensional subscript" i.e. a list with one subscript per
### dimension in the original array.
get_array_block_subscripts <- function(blocks, i, expand.RangeNSBS=FALSE)
{
    nblock <- length(blocks)
    stopifnot(isSingleInteger(i), i >= 1L, i <= nblock)

    ndim <- length(blocks@dim)
    subscripts <- rep.int(alist(foo=), ndim)

    if (blocks@N > ndim)
        return(subscripts)

    i <- i - 1L
    if (blocks@N < ndim) {
        inner_len <- .get_ArrayBlocks_inner_length(blocks)
        i1 <- i %% inner_len
        i2 <- i %/% inner_len
    } else {
        i1 <- i
    }

    k1 <- i1 * blocks@by
    k2 <- k1 + blocks@by
    k1 <- k1 + 1L
    upper_bound <- blocks@dim[[blocks@N]]
    if (k2 > upper_bound)
        k2 <- upper_bound
    if (expand.RangeNSBS) {
        subscript <- k1:k2  # same as doing as.integer() on the RangeNSBS
                            # object below
    } else {
        subscript <- new2("RangeNSBS", subscript=c(k1, k2),
                                       upper_bound=upper_bound,
                                       check=FALSE)
    }
    subscripts[[blocks@N]] <- subscript

    if (blocks@N < ndim) {
        outer_dim <- blocks@dim[(blocks@N + 1L):ndim]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        subscripts[(blocks@N + 1L):ndim] <- as.list(subindex)
    }
    subscripts
}

### NOT exported but used in HDF5Array!
extract_array_block <- function(x, blocks, i)
{
    subscripts <- get_array_block_subscripts(blocks, i,
                                             expand.RangeNSBS=is.array(x))
    subset_by_subscripts(x, subscripts)
}

### NOT exported but used in unit tests.
split_array_in_blocks <- function(x, max_block_len)
{
    blocks <- ArrayBlocks(dim(x), max_block_len)
    lapply(seq_along(blocks),
           function(i) extract_array_block(x, blocks, i))
}

### NOT exported but used in unit tests.
### Rebuild the original array from the subarrays obtained by
### split_array_in_blocks() as an *ordinary* array.
### So if 'x' is an ordinary array, then:
###
###   unsplit_array_from_blocks(split_array_in_blocks(x, max_block_len), x)
###
### should be a no-op for any 'max_block_len' < 'length(x)'.
unsplit_array_from_blocks <- function(subarrays, x)
{
    ans <- combine_array_objects(subarrays)
    dim(ans) <- dim(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks
###
### 3 utility functions to process array-like objects by block.
###

block_APPLY <- function(...)
{
    require_HDF5Array()
    HDF5Array:::block_APPLY(...)
}

block_MAPPLY <- function(...)
{
    require_HDF5Array()
    HDF5Array:::block_MAPPLY(...)
}

block_REDUCE_and_COMBINE <- function(...)
{
    require_HDF5Array()
    HDF5Array:::block_REDUCE_and_COMBINE(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks of columns
###
### 2 convenience wrappers around block_APPLY() and block_REDUCE_and_COMBINE()
### to process a matrix-like object by block of columns.
###

colblock_APPLY <- function(...)
{
    require_HDF5Array()
    HDF5Array:::colblock_APPLY(...)
}

colblock_REDUCE_and_COMBINE <- function(...)
{
    require_HDF5Array()
    HDF5Array:::colblock_REDUCE_and_COMBINE(...)
}

