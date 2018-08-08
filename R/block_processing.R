### =========================================================================
### Block processing utilities
### -------------------------------------------------------------------------
###


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
### read_block() and write_block() generics
###

### Must return an ordinay array and propagate the dimnames.
setGeneric("read_block", signature="x",
    function(x, viewport)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)))
        ans <- standardGeneric("read_block")
        check_returned_array(ans, dim(viewport), "read_block", class(x))
    }
)

### Must return 'x' (possibly modified if it's an in-memory object).
setGeneric("write_block", signature="x",
    function(x, viewport, block)
    {
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), dim(x)),
                  is.array(block),
                  identical(dim(block), dim(viewport)))
        standardGeneric("write_block")
    }
)

### Work on any object 'x' that supports extract_array() e.g. an ordinary
### array or a DelayedArray, HDF5ArraySeed, or DelayedOp object, etc...
### Propagate the dimnames.
setMethod("read_block", "ANY",
    function(x, viewport)
    {
        ## We use extract_array_by_Nindex() instead of extract_array() to
        ## propagate the dimnames.
        Nindex <- makeNindexFromArrayViewport(viewport)
        extract_array_by_Nindex(x, Nindex)
    }
)

setMethod("write_block", "ANY",
    function(x, viewport, block)
    {
        Nindex <- makeNindexFromArrayViewport(viewport)
        replace_by_Nindex(x, Nindex, block)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks
###
### 2 utility functions to process array-like objects by block.
###

.normarg_grid <- function(grid, x)
{
    if (is.null(grid))
        return(blockGrid(x))
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
### of 'x'. If not supplied, the grid returned by 'blockGrid(x)' is used.
### The effective grid (i.e. 'grid' or 'blockGrid(x)'), current block number,
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
            block <- read_block(x, viewport)
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
        block <- read_block(x, viewport)
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

### 'x' is **trusted** to be a logical array-like object.
block_which <- function(x, arr.ind=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(arr.ind))
        stop("'arr.ind' must be TRUE or FALSE")
    ## Return a numeric matrix like one returned by base::arrayInd(), that
    ## is, a matrix where each row is an n-uplet representing an array index.
    FUN <- function(block) {
        b <- currentBlockId(block)
        m <- base::which(block)
        mapToRef(rep.int(b, length(m)), m, effectiveGrid(block), linear=TRUE)
    }
    block_results <- blockApply(x, FUN, grid=grid)
    aind <- do.call(rbind, block_results)
    aind_as_list <- lapply(ncol(aind):1, function(j) aind[ , j])
    oo <- do.call(order, aind_as_list)
    ans <- aind[oo, , drop=FALSE]
    if (!arr.ind)
        ans <- linearInd(ans, dim(x))
    ans
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
                block <- read_block(x, grid[[b]])
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
        ## a trick to make sure that blockGrid() returns linear blocks.
        chunk_grid <- RegularArrayGrid(x_dim, rep.int(1L, length(x_dim)))
    }
    grid <- blockGrid(x, block_maxlen, chunk_grid,
                         block.shape="first-dim-grows-first")
    nblock <- length(grid)
    lapply(seq_len(nblock),
        function(b) {
            if (get_verbose_block_processing())
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            viewport <- grid[[b]]
            block <- read_block(x, viewport)
            block_ans <- APPLY(block, ...)
            if (!is.null(sink)) {
                write_block(sink, viewport, block_ans)
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
    block_maxlen <- max(getDefaultBlockLength(type(x)), x_dim[[1L]])
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

