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
    ## Calling 'defaultSinkAutoGrid()' on 'sink' will produce such grid.
    ## TODO: Note that it might be beneficial to use a grid that is **also**
    ## compatible with the physical chunks of 'x' so we might want to add
    ## that kind of capability to 'defaultSinkAutoGrid()'.
    grid <- defaultSinkAutoGrid(sink)

    FUN <- function(sink, viewport, x, verbose, verbose_read_block)
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
    sinkApply(sink, FUN, x, verbose, verbose_read_block,
              grid=grid, verbose=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### realize()
###

setGeneric("realize", function(x, ...) standardGeneric("realize"))

### Realize an array-like object in memory if 'BACKEND' is NULL, otherwise
### on disk.
### Note that, when 'BACKEND' is not NULL, 'x' gets realized as a "pristine"
### DelayedArray object (e.g. an HDF5Array object), that is, as a DelayedArray
### object that carries no delayed operations. This means that if 'x' is a
### DelayedArray object, the returned object is another DelayedArray object
### that is semantically equivalent to 'x' but where the delayed operations
### carried by 'x' have been realized.
### Will raise an error if 'x' is not an array-like object. Non array-like
### objects (e.g. SummarizedExperiment objects) need to have their own
### realize() method.
setMethod("realize", "ANY",
    function(x, BACKEND=getAutoRealizationBackend())
    {
        if (is.null(dim(x)))
            stop(wmsg("realization of ", class(x)[[1L]], " objects ",
                      "is not supported"))
        if (!is.null(BACKEND)) {
            load_BACKEND_package(BACKEND)
            return(as(x, BACKEND))
        }
        if (is_sparse(x)) {
            as(x, "SparseArraySeed")
        } else {
            ## When 'BACKEND' is NULL, realize() will act as a no-op on an
            ## ordinary matrix or array.
            as.array(x)
        }
    }
)

