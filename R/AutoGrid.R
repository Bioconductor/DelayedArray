### =========================================================================
### Automatic grids
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### defaultAutoGrid()
###

### Guaranteed to return an integer >= 1.
.normarg_block.length <- function(block.length, type)
{
    if (is.null(block.length))
        return(getAutoBlockLength(type))
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
        stop(wmsg("'chunk.grid' must be an ArrayGrid object or NULL"))
    if (!identical(refdim(chunk.grid), dim(x)))
        stop(wmsg("'chunk.grid' is incompatible with 'x'"))
    chunk.grid
}

.normarg_block.shape <- function(block.shape)
{
    if (is.null(block.shape))
        return(getAutoBlockShape())
    if (!(isSingleString(block.shape) &&
          block.shape %in% SUPPORTED_BLOCK_SHAPES))
    {
        in1string <- paste(paste0("\"", SUPPORTED_BLOCK_SHAPES, "\""),
                           collapse=", ")
        stop(wmsg("'block.shape' must be one of ", in1string, ", or NULL"))
    }
    block.shape
}

.block_grid <- function(x_dim, block_len, chunk_grid, block_shape)
{
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
defaultAutoGrid <-
    function(x, block.length=NULL, chunk.grid=NULL, block.shape=NULL)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    block_len <- .normarg_block.length(block.length, type(x))
    chunk_grid <- .normarg_chunk.grid(chunk.grid, x)
    block_shape <- .normarg_block.shape(block.shape)
    .block_grid(x_dim, block_len, chunk_grid, block_shape)
}

blockGrid <- function(...)
{
    .Deprecated("defaultAutoGrid")
    defaultAutoGrid(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowAutoGrid() and colAutoGrid()
###

.get_auto_nrow <- function(x_dim, block.length, x_type)
{
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    block_len <- .normarg_block.length(block.length, x_type)
    nrow <- block_len %/% x_ncol

    if (is.na(nrow)) # NA if x_ncol=0.
        return(x_nrow)
    if (nrow < 1L)
        nrow <- 1L
    min(x_nrow, nrow)
}

### Return a RegularArrayGrid object describing a grid on matrix-like
### object 'x' where the grid elements are blocks made of full rows of 'x'.
rowAutoGrid <- function(x, nrow=NULL, block.length=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (is.null(nrow)) {
        nrow <- .get_auto_nrow(x_dim, block.length, type(x))
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

rowGrid <- function(...)
{
    .Deprecated("rowAutoGrid")
    rowAutoGrid(...)
}

### Return a RegularArrayGrid object describing a grid on matrix-like
### object 'x' where the grid elements are blocks made of full columns of 'x'.
colAutoGrid <- function(x, ncol=NULL, block.length=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (is.null(ncol)) {
        ncol <- .get_auto_nrow(rev(x_dim), block.length, type(x))
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

colGrid <- function(...)
{
    .Deprecated("colAutoGrid")
    colAutoGrid(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### defaultMultAutoGrids()
###

### Performs the same checks as .check_chunkdim() (see chunkGrid.R).
.normarg_z_chunkdim <- function(z_chunkdim, z_dim,
                                what="z", nrow_what="x", ncol_what="y")
{
    if (!is.numeric(z_chunkdim))
        stop(wmsg("'", what, ".chunkdim' must be NULL or an integer vector"))
    if (!is.integer(z_chunkdim))
        z_chunkdim <- as.integer(z_chunkdim)
    if (length(z_chunkdim) != 2L)
        stop(wmsg("'", what, ".chunkdim' must be of length 2"))
    if (S4Vectors:::anyMissingOrOutside(z_chunkdim, 0L))
        stop(wmsg("'", what, ".chunkdim' should not contain negative ",
                  "or NA values"))
    if (z_chunkdim[[1L]] > z_dim[[1L]])
        stop(wmsg("'", what, ".chunkdim[[1L]]' must be ",
                  "<= 'nrow(", nrow_what, ")'"))
    if (z_chunkdim[[1L]] == 0L && z_dim[[1L]] != 0L)
        stop(wmsg("'", what, ".chunkdim[[1L]]' should not be 0 ",
                  "unless 'nrow(", nrow_what, ")' is also 0"))
    if (z_chunkdim[[2L]] > z_dim[[2L]])
        stop(wmsg("'", what, ".chunkdim[[2L]]' must be ",
                  "<= 'ncol(", ncol_what, ")'"))
    if (z_chunkdim[[2L]] == 0L && z_dim[[2L]] != 0L)
        stop(wmsg("'", what, ".chunkdim[[2L]]' should not be 0 ",
                  "unless 'ncol(", ncol_what, ")' is also 0"))
    if (prod(z_chunkdim) > .Machine$integer.max)
        stop(wmsg("The chunk dimensions specified in 'z.chunkdim' are too ",
                  "big. Their product must be <= .Machine$integer.max"))
    z_chunkdim
}

### Define grids for "block matrix multiplication".
###
### The problem
### -----------
### 'x' and 'y' being matrix-like objects that can be multiplied (i.e.
### ncol(x) == nrow(y)), and 'z' being their product (i.e. z = x %*% y),
### defaultMultAutoGrids() tries to solve the problem of finding 3 grids
### (on 'x', 'y', and 'z', respectively) that can be used for "block matrix
### multiplication" of 'x' by 'y'.
### The grid triplet ('x_grid', 'y_grid', 'z_grid') must satisfy:
###   - get_spacings_along(x_grid, 1L)' and 'get_spacings_along(z_grid, 1L)'
###   - get_spacings_along(y_grid, 2L)' and 'get_spacings_along(z_grid, 2L)'
###   - get_spacings_along(x_grid, 2L)' and 'get_spacings_along(y_grid, 1L)'
###
### Many grid triplets satisfy these constraints, from the triplet made of
### the coarsest grids (i.e. each grid is made of a single block covering
### the full matrix e.g. 'x_grid <- RegularArrayGrid(dim(x))') to the triplet
### made of the finest grids (i.e. each grid is made of blocks that contain a
### single matrix element e.g. 'x_grid <- RegularArrayGrid(dim(x), c(1, 1))').
###
### Between these 2 extremes, which would both lead to bad performance in
### general, is a wide range solutions that cover the entire spectrum of
### coarseness. Here is a notable solution with intermediate coarseness:
###   - x_grid <- rowAutoGrid(x, nrow=1)  # blocks on 'x' are single rows
###   - y_grid <- colAutoGrid(y, ncol=1)  # blocks on 'y' are single columns
###   - z_grid <- RegularArrayGrid(dim(z), c(1, 1))
###
### We could also let rowAutoGrid() and colAutoGrid() automatically choose
### the maximum number of rows and columns to put in a block based on the
### "automatic block size" (see ?getAutoBlockSize) by not specifying the
### 'nrow' or 'ncol' arguments:
###  - x_grid <- rowAutoGrid(x)
###  - y_grid <- colAutoGrid(y)
###  - z_grid <- RegularArrayGrid(dim(z), c(a, c))
###    where 'a' is 'nrow(x_grid[[1L]])' and 'c' is 'ncol(y_grid[[1L]])'
###
### Grids that will perform "well"
### ------------------------------
### Amongst all the possible grid triplets, we want to pick up one that will
### perform "well". More precisely:
###
###   (1) In order to cap memory usage during multiplication of any 2 blocks
###       (e.g. 'x_grid[[i, k]]' and 'y_grid[[k, j]]'), we decide (somewhat
###       arbitrarily) that their cumulated length must be as close as possibe
###       to (but must not exceed) half of 'block.length'. The length of the
###       block corresponding to the product ('z_grid[[i, j]]') must also be
###       as close as possibe (but must not exceed) half of 'block.length'.
###
###   (2) We want to minimize the number of times a given block from 'x'
###       and/or 'y' will be loaded in memory. Note that the total number
###       of block multiplications will be 'K * length(z_grid)' where 'K'
###       is 'ncol(x_grid)' (or 'nrow(y_grid)', which is the same). Each
###       block in 'x' will participate to 'ncol(z_grid)' multiplications
###       and each block in 'y' will participate to 'nrow(z_grid)'
###       multiplications. So assuming that blocks are loaded each time
###       they need to participate to a multiplication and not kept in
###       memory, the data in 'x' will be loaded 'ncol(z_grid)' times and
###       the data in 'y' will be loaded 'nrow(z_grid)' times. This means
###       that we want 'z_grid' to be as coarse as possible. Also,
###       depending on which of 'x' or 'y' is bigger, we should try to
###       reduce 'ncol(z_grid)' (i.e. increase coarseness along 'z' 2nd
###       dimension) which would be to the detriment of 'nrow(z_grid)', or
###       vice-versa.
###
###   (3) Ideally we'd want the 3 grids to be compatible with the underlying
###       chunk geometry of the matrix-like object that they are on.
###       Unfortunately this is in general impossible to achieve for the
###       3 grids at the same time (the 3 coarsest grids would achieve this
###       but to the detriment of memory usage) so we'll have to compromise.
###       One thing we don't want to compromise  with though is compatibility
###       of 'z_grid' with the underlying chunk geometry of 'z'. The reason
###       this is that in practice 'z' will typically be a RealizationSink
###       object where the output blocks will be written with write_block()
###       and some realization backends like RleRealizationSink or
###       TENxRealizationSink require write_block() to respect the chunk
###       geometry.
###
### Simplifying the problem
### -----------------------
### In order to simplify the problem we restrict the search space to
### **regular** grids (RegularArrayGrid objects), so no arbitrary grids
### (ArbitraryArrayGrid) will be considered. This reduces the opportunity to
### find the best possible grids but that's ok for now.
### Thanks to this simplification, our problem is reduced to finding 3
### non-negative integer numbers a, b, c from which the the 3 grids will
### be defined as follow:
###  - x_grid <- RegularArrayGrid(dim(x), c(a, b))
###  - y_grid <- RegularArrayGrid(dim(x), c(b, c))
###  - z_grid <- RegularArrayGrid(dim(z), c(a, c))
###

.find_abc <- function(A, B, C, half_block_len,
                      x_chunkdim, y_chunkdim, z_chunkdim)
{
    is_a_multiple_of <- function(m, n) { m == 0L || m %% n == 0L }
    a0 <- z_chunkdim[[1L]]
    if (!is_a_multiple_of(a0, x_chunkdim[[1L]]))
        warning(wmsg("The chunk geometry for 'z' is incompatible with ",
                     "the chunk geometry for 'x' (the number of rows of ",
                     "a chunk in 'z' should be a multiple of the number ",
                     "of rows of a chunk in 'x'). ",
                     "This will likely degrade performance of the block ",
                     "matrix multiplication."), immediate.=TRUE)
    c0 <- z_chunkdim[[2L]]
    if (!is_a_multiple_of(c0, y_chunkdim[[2L]]))
        warning(wmsg("The chunk geometry for 'z' is incompatible with ",
                     "the chunk geometry for 'y' (the number of columns of ",
                     "a chunk in 'z' should be a multiple of the number ",
                     "of columns of a chunk in 'y'). ",
                     "This will likely degrade performance of the block ",
                     "matrix multiplication."), immediate.=TRUE)
    if (C > A && A != 0L)
        a0 <- a0 * as.integer(C / A)
    if (A > C && C != 0L)
        c0 <- c0 * as.integer(A / C)

    ## Find 'a' and 'c'.
    z_dim <- c(A, C)
    z_tile_dim <- pmin(c(a0, c0), z_dim)
    z_tile_grid <- RegularArrayGrid(z_dim, z_tile_dim)
    z_grid <- .block_grid(z_dim, max(as.integer(half_block_len), 1L),
                          z_tile_grid, "scale")
    z_spacings <- dim(z_grid[[1L]])
    a <- z_spacings[[1L]]
    c <- z_spacings[[2L]]

    ## Find 'b'.
    ## TODO: Revisit this.
    b <- max(as.integer(half_block_len / (a + c)), 1L)
    b0 <- max(x_chunkdim[[2L]], y_chunkdim[[1L]])
    r <- b %% b0
    if (r != 0L)
        b <- b + (b0 - r)  # increase 'b' to the next multiple of 'b0'
    if (b > B)
        b <- B

    c(a=a, b=b, c=c)
}

### Return list(x_grid, y_grid, z_grid).
defaultMultAutoGrids <-
    function(x, y, block.length=NULL,
             x.chunkdim=NULL, y.chunkdim=NULL, z.chunkdim=NULL)
{
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (length(x_dim) != 2L || length(y_dim) != 2L)
        stop(wmsg("'x' and 'y' must have exactly 2 dimensions"))
    x_ncol <- x_dim[[2L]]
    y_nrow <- y_dim[[1L]]
    if (x_ncol != y_nrow)
        stop(wmsg("'ncol(x)' must be equal to 'nrow(y)'"))
    A <- x_dim[[1L]]
    B <- x_ncol  # same as 'y_nrow'
    C <- y_dim[[2L]]
    z_dim <- c(A, C)

    ## Assuming that 'x', 'y', and 'z' are all of type "double". This could
    ## be refined.
    half_block_len <- .normarg_block.length(block.length, "double") / 2

    if (is.null(x.chunkdim)) {
        x_chunkdim <- chunkdim(x)
        if (is.null(x_chunkdim))
            x_chunkdim <- pmin(1L, x_dim)
    } else {
        x_chunkdim <- .normarg_z_chunkdim(x.chunkdim, x_dim,
                                          what="x", ncol_what="x")
    }
    if (is.null(y.chunkdim)) {
        y_chunkdim <- chunkdim(y)
        if (is.null(y_chunkdim))
            y_chunkdim <- pmin(1L, y_dim)
    } else {
        y_chunkdim <- .normarg_z_chunkdim(y.chunkdim, y_dim,
                                          what="y", nrow_what="y")
    }
    if (is.null(z.chunkdim)) {
        z_chunkdim <- pmin(1L, z_dim)
    } else {
        z_chunkdim <- .normarg_z_chunkdim(z.chunkdim, z_dim)
    }

    abc <- .find_abc(A, B, C, half_block_len,
                     x_chunkdim, y_chunkdim, z_chunkdim)

    x_grid <- RegularArrayGrid(x_dim, abc[-3L])
    y_grid <- RegularArrayGrid(y_dim, abc[-1L])
    z_grid <- RegularArrayGrid(z_dim, abc[-2L])
    list(x_grid, y_grid, z_grid)
}

multGrids <- function(...)
{
    .Deprecated("defaultMultAutoGrids")
    defaultMultAutoGrids(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoGridMaker()
###

### We set the automatic grid maker to defaultAutoGrid() by default.
setAutoGridMaker <- function(GRIDMAKER="defaultAutoGrid")
{
    match.fun(GRIDMAKER)  # sanity check
    set_user_option("auto.grid.maker", GRIDMAKER)
}

getAutoGridMaker <- function() get_user_option("auto.grid.maker")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_grid()
###

normarg_grid <- function(grid, x)
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

