### =========================================================================
### realize()
### -------------------------------------------------------------------------
###


setGeneric("realize", function(x) standardGeneric("realize"))

setMethod("realize", "ANY",
    function(x) as(DelayedArray(x), getRealizeBackend())
)

