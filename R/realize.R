### =========================================================================
### realize()
### -------------------------------------------------------------------------
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

