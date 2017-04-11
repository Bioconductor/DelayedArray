### =========================================================================
### ArrayBlocks objects
### -------------------------------------------------------------------------
###


setClass("ArrayBlocks",
    contains="List",
    representation(
        dim="integer",
        max_block_len="integer",
        N="integer",
        by="integer"
    ),
    prototype(elementType="list")
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
    new("ArrayBlocks", dim=dim, max_block_len=max_block_len, N=N, by=by)
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

### Return an IRanges object with 1 range per dimension.
get_block_ranges <- function(blocks, i)
{
    nblock <- length(blocks)
    stopifnot(isSingleInteger(i), i >= 1L, i <= nblock)

    ndim <- length(blocks@dim)
    ans <- IRanges(rep.int(1L, ndim), blocks@dim)

    if (blocks@N > ndim)
        return(ans)

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
    start(ans)[[blocks@N]] <- k1
    end(ans)[[blocks@N]] <- k2

    if (blocks@N < ndim) {
        outer_dim <- blocks@dim[(blocks@N + 1L):ndim]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        ans[(blocks@N + 1L):ndim] <- IRanges(subindex, width=1L)
    }
    ans
}

get_array_block_subscripts <- function(blocks, i)
{
    make_subscripts_from_block_ranges(get_block_ranges(blocks, i), blocks@dim)
}

setMethod("getListElement", "ArrayBlocks",
    function(x, i, exact=TRUE)
    {
        i <- normalizeDoubleBracketSubscript(i, x, exact=exact, 
                                             error.if.nomatch=TRUE)
        get_array_block_subscripts(x, i)
    }
)

setMethod("show", "ArrayBlocks",
    function(object)
    {
        dim_in1string <- paste0(object@dim, collapse=" x ")
        cat(class(object), " object with ", length(object), " blocks ",
            "of length <= ", object@max_block_len, " on a ",
            dim_in1string, " array:\n", sep="")
        for (i in seq_along(object)) {
            subscripts <- object[[i]]
            cat("[[", i, "]]: [", subscripts_as_string(subscripts), "]\n",
                sep="")
        }
    }
)

extract_array_block <- function(x, blocks, i)
{
    subscripts <- get_array_block_subscripts(blocks, i)
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

