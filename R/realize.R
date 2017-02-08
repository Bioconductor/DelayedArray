### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


setGeneric("realize", function(x) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x)
    {
        ans <- as(DelayedArray(x), getRealizeBackend())
        ## Temporarily needed because coercion to HDF5Array currently drops
        ## the dimnames. See R/writeHDF5Array.R in the HDF5Array package for
        ## more information about this.
        ## TODO: Remove line below when this is addressed.
        dimnames(ans) <- dimnames(x)
        ans
    }
)

