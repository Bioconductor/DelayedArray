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

### Return an ArrayBlocks object i.e. a collection of viewports on a reference
### array. This collection as the following properties:
###   (a) The collection of viewports defines a partitioning of the reference
###       array i.e. they fully cover it and don't overlap each other.
###   (b) Each viewport delimits a block that has a length (i.e. nb of
###       elements) <= 'max_block_len'.
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

.get_ArrayBlocks_inner_length <- function(grid)
{
    ndim <- length(grid@dim)
    if (grid@N > ndim)
        return(if (any(grid@dim == 0L)) 0L else 1L)
    inner_len <- grid@dim[[grid@N]] %/% grid@by
    by2 <- grid@dim[[grid@N]] %% grid@by
    if (by2 != 0L)  # 'grid' contains truncated viewports
        inner_len <- inner_len + 1L
    inner_len
}

.get_ArrayBlocks_outer_length <- function(grid)
{
    ndim <- length(grid@dim)
    if (grid@N >= ndim)
        return(1L)
    outer_dim <- grid@dim[(grid@N + 1L):ndim]
    prod(outer_dim)
}

### Return the number of viewports in grid 'x'.
setMethod("length", "ArrayBlocks",
    function(x)
        .get_ArrayBlocks_inner_length(x) * .get_ArrayBlocks_outer_length(x)
)

get_block_lengths <- function(grid)
{
    p <- prod(grid@dim[seq_len(grid@N - 1L)])
    ndim <- length(grid@dim)
    if (grid@N > ndim)
        return(p)
    fb_len <- p * grid@by  # full block length
    lens <- rep.int(fb_len, grid@dim[[grid@N]] %/% grid@by)
    by2 <- grid@dim[[grid@N]] %% grid@by
    if (by2 != 0L) {         # 'grid' contains truncated viewports
        tb_len <- p * by2    # truncated block length
        lens <- c(lens, tb_len)
    }
    rep.int(lens, .get_ArrayBlocks_outer_length(grid))
}

### Return an ArrayBlock object.
setMethod("getListElement", "ArrayBlocks",
    function(x, i, exact=TRUE)
    {
        i <- normalizeDoubleBracketSubscript(i, x, exact=exact, 
                                             error.if.nomatch=TRUE)

        ## We start with a block that covers the whole array.
        ans <- ArrayViewport(x@dim)

        ndim <- length(x@dim)
        if (x@N > ndim)
            return(ans)

        i <- i - 1L
        if (x@N < ndim) {
            inner_len <- .get_ArrayBlocks_inner_length(x)
            i1 <- i %% inner_len
            i2 <- i %/% inner_len
        } else {
            i1 <- i
        }

        k1 <- i1 * x@by
        k2 <- k1 + x@by
        k1 <- k1 + 1L
        upper_bound <- x@dim[[x@N]]
        if (k2 > upper_bound)
            k2 <- upper_bound
        start(ans@ranges)[[x@N]] <- k1
        end(ans@ranges)[[x@N]] <- k2

        if (x@N < ndim) {
            outer_dim <- x@dim[(x@N + 1L):ndim]
            subindex <- arrayInd(i2 + 1L, outer_dim)
            ans@ranges[(x@N + 1L):ndim] <- IRanges(subindex, width=1L)
        }
        ans
    }
)

setMethod("show", "ArrayBlocks",
    function(object)
    {
        dim_in1string <- paste0(object@dim, collapse=" x ")
        cat(class(object), " object with ", length(object), " viewports ",
            "of length <= ", object@max_block_len, " on a ",
            dim_in1string, " array:\n", sep="")
        for (i in seq_along(object)) {
            s <- make_string_from_ArrayViewport(object[[i]], with.brackets=TRUE)
            cat("[[", i, "]]: ", s, "\n", sep="")
        }
    }
)

extract_array_block <- function(x, grid, i)
{
    Nindex <- makeNindexFromArrayViewport(grid[[i]])
    subset_by_Nindex(x, Nindex)
}

### NOT exported but used in unit tests.
split_array_in_blocks <- function(x, max_block_len)
{
    grid <- ArrayBlocks(dim(x), max_block_len)
    lapply(seq_along(grid),
           function(i) extract_array_block(x, grid, i))
}

### NOT exported but used in unit tests.
### Rebuild the original array from the blocks obtained by
### split_array_in_blocks() as an *ordinary* array.
### So if 'x' is an ordinary array, then:
###
###   unsplit_array_from_blocks(split_array_in_blocks(x, max_block_len), x)
###
### should be a no-op for any 'max_block_len' < 'length(x)'.
unsplit_array_from_blocks <- function(blocks, x)
{
    ans <- combine_array_objects(blocks)
    dim(ans) <- dim(x)
    ans
}

