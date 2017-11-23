### =========================================================================
### ConformableSeedCombiner objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("ConformableSeedCombiner",
    representation(
        seeds="list",             # List of n conformable array-like objects
                                  # to combine. Each object is expected to
                                  # satisfy the "seed contract" i.e. to
                                  # support dim(), dimnames(), and
                                  # extract_array().

        COMBINING_OP="function",  # n-ary operator to combine the seeds.

        Rargs="list"              # Additional arguments to the n-ary
                                  # operator.
    ),
    prototype(
        seeds=list(new("array")),
        COMBINING_OP=identity
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

.validate_ConformableSeedCombiner <- function(x)
{
    ## 'seeds' slot.
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))
    if (!.objects_are_conformable_arrays(x@seeds))
        return(wmsg2("'x@seeds' must be a list of conformable ",
                     "array-like objects"))
    TRUE
}

setValidity2("ConformableSeedCombiner", .validate_ConformableSeedCombiner)

new_ConformableSeedCombiner <- function(seed=new("array"), ...,
                                        COMBINING_OP=identity,
                                        Rargs=list())
{
    seeds <- unname(list(seed, ...))
    COMBINING_OP <- match.fun(COMBINING_OP)
    new2("ConformableSeedCombiner", seeds=seeds,
                                    COMBINING_OP=COMBINING_OP,
                                    Rargs=Rargs)
}

### Implement the "seed contract" i.e. dim(), dimnames(), and
### extract_array().

.get_ConformableSeedCombiner_dim <- function(x) dim(x@seeds[[1L]])

setMethod("dim", "ConformableSeedCombiner",
    .get_ConformableSeedCombiner_dim
)

.get_ConformableSeedCombiner_dimnames <- function(x)
{
    combine_dimnames(x@seeds)
}

setMethod("dimnames", "ConformableSeedCombiner",
    .get_ConformableSeedCombiner_dimnames
)

.extract_array_from_ConformableSeedCombiner <- function(x, index)
{
    arrays <- lapply(x@seeds, extract_array, index)
    do.call(x@COMBINING_OP, c(arrays, x@Rargs))
}

setMethod("extract_array", "ConformableSeedCombiner",
    .extract_array_from_ConformableSeedCombiner
)

