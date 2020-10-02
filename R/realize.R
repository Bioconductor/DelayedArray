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

BLOCK_write_to_sink <- function(sink, x, verbose=NA)
{
    stopifnot(identical(dim(x), dim(sink)))
    verbose <- normarg_verbose(verbose)

    ## 'x' and 'sink' might both have their physical chunks but we must
    ## choose a grid that is compatible with the physical chunks of 'sink'.
    ## Calling 'defaultAutoGrid()' on 'sink' will produce such grid.
    ## Note that it might be beneficial to use a grid that is also compatible
    ## with the physical chunks of 'x' so we might want to come up with a
    ## dedicated utility for that e.g. 'defaultAutoGrid2(sink, x)'.
    ## Also by using block.shape="first-dim-grows-first" in the call below
    ## we'll get a grid that guarentees linear writing to the sink in case
    ## 'chunkdim(sink)' is NULL.
    grid <- defaultAutoGrid(sink, block.shape="first-dim-grows-first")

    FUN <- function(viewport, sink, x, verbose, verbose_read_block)
    {
        effective_grid <- effectiveGrid()
        current_block_id <- currentBlockId()
        if (verbose) {
            x_is_sparse <- is_sparse(x)
            nblock <- length(effective_grid)
            block <- verbose_read_block(x, viewport, x_is_sparse,
                                        as.sparse=NA, current_block_id, nblock)
        } else {
            block <- read_block(x, viewport, as.sparse=NA)
        }
        if (verbose)
            message("\\ Writing it ... ", appendLF=FALSE)
        sink <- write_block(sink, viewport, block)
        if (verbose)
            message("OK")
        sink
    }
    viewportReduce(FUN, grid, sink,
                   x, verbose, verbose_read_block,
                   verbose=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

setGeneric("realize", function(x, ...) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x, BACKEND=getAutoRealizationBackend())
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

