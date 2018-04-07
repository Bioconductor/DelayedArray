### =========================================================================
### Simplify a tree of DelayedOp objects
### -------------------------------------------------------------------------
###

setGeneric("simplify", function(x) standardGeneric("simplify"))

setMethod("simplify", "ANY", identity)

setMethod("simplify", "DelayedSubset",
    function(x)
    {
        seed <- x@seed
        if (is(seed, "DelayedSubset")) {
            index <- subset_index(seed@index, x@index)
            return(new2("DelayedSubset", seed=seed@seed, index=index))
        }
        if (is(seed, "DelayedAperm")) {
        }
        if (is(seed, "DelayedUnaryIsoOp")) {
        }
        if (is(seed, "DelayedDimnames")) {
        }
        x
    }
)

