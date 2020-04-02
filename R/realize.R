### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Internal helper to support block by block realization
###
### Used by coercion to RleArray, HDF5Array::writeHDF5Array(), and
### HDF5Array::writeTENxMatrix().
###

BLOCK_write_to_sink <- function(x, sink)
{
    stopifnot(identical(dim(x), dim(sink)))
    ## 'x' and 'sink' might both have their physical chunks but we must
    ## choose a grid that is compatible with the physical chunks of 'sink'.
    ## Calling 'blockGrid()' on 'sink' will produce such grid.
    ## Note that it might be beneficial to use a grid that is also compatible
    ## with the physical chunks of 'x' so we might want to come up with a
    ## dedicated utility for that e.g. 'blockGrid2(sink, x)'.
    ## Also by using block.shape="first-dim-grows-first" in the call below
    ## we'll get a grid that guarentees linear writing to the sink in case
    ## 'chunkdim(sink)' is NULL.
    grid <- blockGrid(sink, block.shape="first-dim-grows-first")
    nblock <- length(grid)
    x_is_sparse <- is_sparse(x)
    for (bid in seq_len(nblock)) {
        viewport <- grid[[bid]]
        if (x_is_sparse) {
            if (get_verbose_block_processing())
                message("Realizing sparse block ", bid, "/", nblock, " ... ",
                        appendLF=FALSE)
            sparse_block <- read_sparse_block(x, viewport)
            if (get_verbose_block_processing())
                message("OK, writing it ... ",
                        appendLF=FALSE)
            write_sparse_block(sink, viewport, sparse_block)
            if (get_verbose_block_processing())
                message("OK")
        } else {
            if (get_verbose_block_processing())
                message("Realizing block ", bid, "/", nblock, " ... ",
                        appendLF=FALSE)
            block <- read_block(x, viewport)
            if (get_verbose_block_processing())
                message("OK, writing it ... ",
                        appendLF=FALSE)
            write_block(sink, viewport, block)
            if (get_verbose_block_processing())
                message("OK")
        }
    }
    sink
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

setGeneric("realize", function(x, ...) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x, BACKEND=getRealizationBackend())
    {
        x <- DelayedArray(x)
        if (is.null(BACKEND))
            return(DelayedArray(as.array(x)))
        load_BACKEND_package(BACKEND)
        ans <- as(x, BACKEND)
        ## Temporarily needed because coercion to HDF5Array currently drops
        ## the dimnames. See R/writeHDF5Array.R in the HDF5Array package for
        ## more information about this.
        ## TODO: Remove line below when this is addressed.
        set_dimnames(ans, dimnames(x))
    }
)

