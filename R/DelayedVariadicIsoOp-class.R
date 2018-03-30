### =========================================================================
### DelayedVariadicIsoOp objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("DelayedVariadicIsoOp",
    representation(
        seeds="list",   # List of conformable array-like objects to combine.
                        # Each object is expected to satisfy the "seed
                        # contract" i.e. to support dim(), dimnames(), and
                        # extract_array().

        OP="function",  # The function to use to combine the seeds. It should
                        # act as an isomorphism i.e. always return an array
                        # parallel to the input arrays (i.e. same dimensions).

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        seeds=list(new("array")),
        OP=identity
    )
)

.objects_are_conformable_arrays <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    first_ndim <- ndims[[1L]]
    if (!all(ndims == first_ndim))
        return(FALSE)
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return(FALSE)
    dims <- matrix(tmp, nrow=first_ndim)
    first_dim <- dims[ , 1L]
    all(dims == first_dim)
}

.validate_DelayedVariadicIsoOp <- function(x)
{
    ## 'seeds' slot.
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))
    if (!.objects_are_conformable_arrays(x@seeds))
        return(wmsg2("'x@seeds' must be a list of conformable ",
                     "array-like objects"))
    TRUE
}

setValidity2("DelayedVariadicIsoOp", .validate_DelayedVariadicIsoOp)

new_DelayedVariadicIsoOp <- function(seed=new("array"), ...,
                                     OP=identity, Rargs=list())
{
    seeds <- unname(list(seed, ...))
    OP <- match.fun(OP)
    new2("DelayedVariadicIsoOp", seeds=seeds, OP=OP, Rargs=Rargs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array()
###

.get_DelayedVariadicIsoOp_dim <- function(x) dim(x@seeds[[1L]])

setMethod("dim", "DelayedVariadicIsoOp",
    .get_DelayedVariadicIsoOp_dim
)

.get_DelayedVariadicIsoOp_dimnames <- function(x) combine_dimnames(x@seeds)

setMethod("dimnames", "DelayedVariadicIsoOp",
    .get_DelayedVariadicIsoOp_dimnames
)

.extract_array_from_DelayedVariadicIsoOp <- function(x, index)
{
    arrays <- lapply(x@seeds, extract_array, index)
    do.call(x@OP, c(arrays, x@Rargs))
}

setMethod("extract_array", "DelayedVariadicIsoOp",
    .extract_array_from_DelayedVariadicIsoOp
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###
### DelayedVariadicIsoOp objects don't support the path() getter/setter for
### now.
###

setMethod("path", "DelayedVariadicIsoOp",
    function(object, ...)
        stop(wmsg("path() is not supported on a DelayedArray ",
                  "object with multiple leaf seeds at the moment"))
)

setReplaceMethod("path", "DelayedVariadicIsoOp",
    function(object, ..., value)
        stop(wmsg("the path() setter is not supported on a DelayedArray ",
                  "object with multiple leaf seeds at the moment"))
)

