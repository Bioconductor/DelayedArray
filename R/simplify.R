### =========================================================================
### Simplify a tree of DelayedOp objects
### -------------------------------------------------------------------------
###

setGeneric("simplify", function(x) standardGeneric("simplify"))

setMethod("simplify", "ANY", identity)

setMethod("simplify", "DelayedSubset",
    function(x)
    {
        x1 <- x@seed
        if (is(x1, "DelayedSubset")) {
            ## SQUASH
            index <- subset_index(x1@index, x@index)
            x <- new2("DelayedSubset", seed=x1@seed, index=index)
            return(x)
        }
        if (is(x1, "DelayedAperm")) {
            ## SWAP
            index2 <- rep.int(list(NULL), length(dim(x1@seed)))
            index2[x1@perm] <- x@index
            x2 <- new2("DelayedSubset", seed=x1@seed, index=index2)
            x1@seed <- simplify(x2)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOp")) {
            ## SWAP
            Largs <- subset_args(x1@Largs, x1@Lalong, x@index)
            Rargs <- subset_args(x1@Rargs, x1@Ralong, x@index)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x,
                                                  Largs=Largs,
                                                  Rargs=Rargs)
            return(x1)
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x_dimnames <- dimnames(x)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- new_DelayedDimnames(x, x_dimnames)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedAperm",
    function(x)
    {
        x1 <- x@seed
        if (is(x1, "DelayedAperm")) {
            ## SQUASH + REMOVE IF NO-OP
            x1@perm <- x1@perm[x@perm]
            if (isNoOp(x1))
                return(x1@seed)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOp")) {
            ## SWAP
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x_dimnames <- dimnames(x)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- new_DelayedDimnames(x, x_dimnames)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedUnaryIsoOp",
    function(x)
    {
        x1 <- x@seed
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x@seed <- x1@seed
            x1@seed <- x
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedDimnames",
    function(x)
    {
        x1 <- x@seed
        if (is(x1, "DelayedDimnames")) {
            ## SQUASH + REMOVE IF NO-OP
            x <- new_DelayedDimnames(x1@seed, dimnames(x))
            if (isNoOp(x))
                return(x@seed)
            return(x)
        }
        x
    }
)

