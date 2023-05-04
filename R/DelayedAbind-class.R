### =========================================================================
### DelayedAbind objects
### -------------------------------------------------------------------------
###
### Representation of a delayed abind() operation.
###

setClass("DelayedAbind",
    contains="DelayedNaryOp",
    representation(
        along="integer"  # Single integer indicating the dimension along
                         # which to bind the input objects.
    ),
    prototype(
        along=1L
    )
)

.validate_DelayedAbind <- function(x)
{
    ## 'along' slot.
    if (!(isSingleInteger(x@along) && x@along >= 1L))
        return("'x@along' must be a single positive integer")
    ndim <- length(dim(x@seeds[[1L]]))
    if (ndim < x@along)
        return(paste0("the array-like objects to bind must have at least ",
                      x@along, " dimensions for this binding operation"))

    dims <- S4Arrays:::get_dims_to_bind(x@seeds, x@along)
    if (is.character(dims))
        return(dims)
    TRUE
}

setValidity2("DelayedAbind", .validate_DelayedAbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

new_DelayedAbind <- function(seeds, along)
{
    new2("DelayedAbind", seeds=seeds, along=along)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_noop() method
###

setMethod("is_noop", "DelayedAbind",
    function(x) length(x@seeds) == 1L
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedAbind

.DelayedAbind_summary <-
    function(object) sprintf("Abind (along=%d)", object@along)

summary.DelayedAbind <-
    function(object, ...) .DelayedAbind_summary(object, ...)

setMethod("summary", "DelayedAbind", summary.DelayedAbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###

.get_DelayedAbind_dim <- function(x)
{
    dims <- S4Arrays:::get_dims_to_bind(x@seeds, x@along)
    S4Arrays:::combine_dims_along(dims, x@along)
}

setMethod("dim", "DelayedAbind", .get_DelayedAbind_dim)

.get_DelayedAbind_dimnames <- function(x)
{
    dims <- S4Arrays:::get_dims_to_bind(x@seeds, x@along)
    S4Arrays:::combine_dimnames_along(x@seeds, dims, x@along)
}

setMethod("dimnames", "DelayedAbind", .get_DelayedAbind_dimnames)

.extract_array_from_DelayedAbind <- function(x, index)
{
    i <- index[[x@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        arrays <- lapply(x@seeds, extract_array, index)
        ## Bind the ordinary arrays in 'arrays'.
        ans <- do.call(S4Arrays:::simple_abind, c(arrays, list(along=x@along)))
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- S4Arrays:::get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- S4Arrays:::get_part_index(i, breakpoints)
    split_part_idx <- S4Arrays:::split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[x@along]] <- split_part_idx[[s]]
        extract_array(x@seeds[[s]], index)
    }
    arrays <- lapply(seq_along(x@seeds), FUN)

    ## Bind the ordinary arrays in 'arrays'.
    ans <- do.call(S4Arrays:::simple_abind, c(arrays, list(along=x@along)))

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector("list", length=length(index))
    Nindex[[x@along]] <- S4Arrays:::get_rev_index(part_idx)
    extract_array(ans, Nindex)
}

setMethod("extract_array", "DelayedAbind", .extract_array_from_DelayedAbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedAbind",
    function(x)
    {
        all(vapply(x@seeds, is_sparse, logical(1), USE.NAMES=FALSE))
    }
)

.OLD_extract_sparse_array_DelayedAbind <- function(x, index)
{
    i <- index[[x@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        arrays <- lapply(x@seeds, OLD_extract_sparse_array, index)
        ## Bind the SparseArraySeed objects in 'arrays'.
        ans <- abind_SparseArraySeed_objects(arrays, x@along)
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- S4Arrays:::get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- S4Arrays:::get_part_index(i, breakpoints)
    split_part_idx <- S4Arrays:::split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[x@along]] <- split_part_idx[[s]]
        OLD_extract_sparse_array(x@seeds[[s]], index)
    }
    arrays <- lapply(seq_along(x@seeds), FUN)

    ## Bind the SparseArraySeed objects in 'arrays'.
    ans <- abind_SparseArraySeed_objects(arrays, x@along)

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector("list", length=length(index))
    Nindex[[x@along]] <- S4Arrays:::get_rev_index(part_idx)
    OLD_extract_sparse_array(ans, Nindex)
}

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "OLD_extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("OLD_extract_sparse_array", "DelayedAbind",
    .OLD_extract_sparse_array_DelayedAbind
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Backward compatibility with DelayedArray < 0.5.24
###
### In DelayedArray 0.5.24 the SeedBinder class got renamed DelayedAbind.
### DelayedArray objects serialized with DelayedArray < 0.5.24 might contain
### SeedBinder instances nested in their "seed" slot so we need to keep
### the class around for now.
###

setClass("SeedBinder", contains="DelayedAbind")

setMethod("updateObject", "SeedBinder",
    function(object, ..., verbose=FALSE)
    {
        class(object) <- "DelayedAbind"
        callNextMethod()
    }
)

