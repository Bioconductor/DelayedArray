### =========================================================================
### DelayedUnaryIsoOp objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("DelayedUnaryIsoOp",
    representation(
        seed="ANY",     # An array-like object expected to satisfy the
                        # "seed contract" i.e. to support dim(), dimnames(),
                        # and extract_array().

        OP="function",  # The function to apply to the seed (e.g. `+` or
                        # log). It should act as an isomorphism i.e. always
                        # return an array parallel to the input array (i.e.
                        # same dimensions).

        Largs="list",   # Additional left arguments to OP.

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        seed=new("array"),
        OP=identity
    )
)

new_DelayedUnaryIsoOp <- function(seed=new("array"),
                                  OP=identity, Largs=list(), Rargs=list())
{
    OP <- match.fun(OP)
    new2("DelayedUnaryIsoOp", seed=seed, OP=OP, Largs=Largs, Rargs=Rargs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array()
###

.get_DelayedUnaryIsoOp_dim <- function(x) dim(x@seed)

setMethod("dim", "DelayedUnaryIsoOp", .get_DelayedUnaryIsoOp_dim)

.get_DelayedUnaryIsoOp_dimnames <- function(x) dimnames(x@seed)

setMethod("dimnames", "DelayedUnaryIsoOp", .get_DelayedUnaryIsoOp_dimnames)

.extract_array_from_DelayedUnaryIsoOp <- function(x, index)
{
    a <- extract_array(x@seed, index)
    do.call(x@OP, c(x@Largs, list(a), x@Rargs))
}

setMethod("extract_array", "DelayedUnaryIsoOp",
    .extract_array_from_DelayedUnaryIsoOp
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seed() getter/setter
###

setMethod("seed", "DelayedUnaryIsoOp", function(x) x@seed)

setReplaceMethod("seed", "DelayedUnaryIsoOp",
    function(x, value)
    {
        x@seed <- normalize_seed_replacement_value(value, seed(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### The path of a DelayedUnaryIsoOp object is the path of its seed so path()
### will work only on a DelayedUnaryIsoOp object with a seed that supports
### path().
### For example it will work if the seed is an on-disk object (e.g. an
### HDF5ArraySeed object) but not if it's an in-memory object (e.g. an
### ordinary array or RleArraySeed object).
setMethod("path", "DelayedUnaryIsoOp",
    function(object, ...) path(seed(object), ...)
)

### The path() setter sets the supplied path on the seed of the
### DelayedUnaryIsoOp object so it will work out-of-the-box on any
### DelayedUnaryIsoOp object with a seed that supports the path() setter.
### For example it will work if the seed is an HDF5ArraySeed object.
setReplaceMethod("path", "DelayedUnaryIsoOp",
    function(object, ..., value)
    {
        path(seed(object), ...) <- value
        object
    }
)

