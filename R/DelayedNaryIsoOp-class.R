### =========================================================================
### DelayedNaryIsoOp objects
### -------------------------------------------------------------------------
###
### Representation of a delayed N-ary isometric operation.
### The input arrays must be "conformable" i.e. they all must have the same
### dimensions.
###

setClass("DelayedNaryIsoOp",
    contains="DelayedNaryOp",
    representation(
        OP="function",  # The function to use to combine the input objects.
                        # Should act as an isomorphism i.e. always return an
                        # array-like object **parallel** to the input objects
                        # (i.e. with the same dimensions).

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        OP=identity
    )
)

.arrays_are_conformable <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    first_ndim <- ndims[[1L]]
    if (!all(ndims == first_ndim))
        return(FALSE)
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return(FALSE)
    dims <- matrix(tmp, ncol=length(objects))
    first_dim <- dims[ , 1L]
    all(dims == first_dim)
}

.validate_DelayedNaryIsoOp <- function(x)
{
    ## 'seeds' slot.
    if (!.arrays_are_conformable(x@seeds))
        return("'x@seeds' must be a list of conformable array-like objects")
    TRUE
}

setValidity2("DelayedNaryIsoOp", .validate_DelayedNaryIsoOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

new_DelayedNaryIsoOp <- function(OP=identity, seed=new("array"), ...,
                                 Rargs=list())
{
    OP <- match.fun(OP)
    seeds <- unname(list(seed, ...))
    if (!.arrays_are_conformable(seeds))
        stop(wmsg("non-conformable array-like objects"))
    new2("DelayedNaryIsoOp", seeds=seeds, OP=OP, Rargs=Rargs, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedNaryIsoOp

.DelayedNaryIsoOp_summary <- function(object) "N-ary iso op"

summary.DelayedNaryIsoOp <-
    function(object, ...) .DelayedNaryIsoOp_summary(object, ...)

setMethod("summary", "DelayedNaryIsoOp", summary.DelayedNaryIsoOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###

setMethod("dim", "DelayedNaryIsoOp", function(x) dim(x@seeds[[1L]]))

setMethod("dimnames", "DelayedNaryIsoOp",
    function(x) S4Arrays:::get_first_non_NULL_dimnames(x@seeds)
)

setMethod("extract_array", "DelayedNaryIsoOp",
    function(x, index)
    {
        arrays <- lapply(x@seeds, extract_array, index)
        do.call(x@OP, c(arrays, x@Rargs))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedNaryIsoOp",
    function(x)
    {
        ok <- vapply(x@seeds, is_sparse, logical(1), USE.NAMES=FALSE)
        if (!all(ok))
            return(FALSE)
        if (length(x@Rargs) != 0L)
            return(FALSE)
        ## Structural sparsity will be propagated if the operation in
        ## x@OP preserves the zeros. To find out whether zeros are preserved
        ## or not, we replace each current seed with an array of one "zero",
        ## that is, with an ordinary array of the same number of dimensions
        ## and type as the seed, but with a single "zero" element. Then we
        ## apply the n-ary operation in x@OP to them and see whether the
        ## zero were preserved or not.
        seed_ndim <- length(dim(x@seeds[[1L]]))
        x@seeds <- lapply(x@seeds,
            function(seed) S4Arrays:::make_one_zero_array(type(seed), seed_ndim))
        ## Same as 'as.array(x)' but doesn't try to propagate the dimnames.
        a0 <- extract_array(x, vector("list", length=seed_ndim))
        S4Arrays:::is_filled_with_zeros(a0)
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "OLD_extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("OLD_extract_sparse_array", "DelayedNaryIsoOp",
    function(x, index)
    {
        stop("NOT IMPLEMENTED YET!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Backward compatibility with DelayedArray < 0.5.24
###
### In DelayedArray 0.5.24 the ConformableSeedCombiner class got renamed
### DelayedNaryIsoOp. DelayedArray objects serialized with DelayedArray <
### 0.5.24 might contain ConformableSeedCombiner instances nested in their
### "seed" slot so we need to keep the class around for now.
###

setClass("ConformableSeedCombiner", contains="DelayedNaryIsoOp")

setMethod("updateObject", "ConformableSeedCombiner",
    function(object, ..., verbose=FALSE)
    {
        object <- new2("DelayedNaryIsoOp", seeds=object@seeds,
                                           OP=object@COMBINING_OP,
                                           Rargs=object@Rargs)
        callNextMethod()
    }
)

