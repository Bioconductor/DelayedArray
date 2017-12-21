### =========================================================================
### SeedDimPicker objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("SeedDimPicker",
    representation(
        seed="ANY",            # An array-like object expected to satisfy
                               # the "seed contract" i.e. to support dim(),
                               # dimnames(), and extract_array().

        picked_dims="integer"  # Index into dim(seed) specifying the seed
                               # dimensions to keep.
    ),
    prototype(
        seed=new("array"),
        picked_dims=1L
    )
)

.validate_SeedDimPicker <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)
    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))
    ## 'picked_dims' slot.
    if (length(x@picked_dims) == 0L)
        return(wmsg2("'x@picked_dims' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@picked_dims, 1L, seed_ndim))
        return(wmsg2("all values in 'x@picked_dims' must be >= 1 ",
                     "and <= 'seed_ndim'"))
    if (anyDuplicated(x@picked_dims))
        return(wmsg2("'x@picked_dims' cannot have duplicates"))
    if (!all(seed_dim[-x@picked_dims] == 1L))
        return(wmsg2("dimensions to drop from 'x' must be equal to 1"))
    TRUE
}

setValidity2("SeedDimPicker", .validate_SeedDimPicker)

new_SeedDimPicker <- function(seed, picked_dims)
{
    seed <- remove_pristine_DelayedArray_wrapping(seed)
    new2("SeedDimPicker", seed=seed, picked_dims=picked_dims)
}

### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array().

.get_SeedDimPicker_dim <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_dim[x@picked_dims]
}

setMethod("dim", "SeedDimPicker", .get_SeedDimPicker_dim)

.get_SeedDimPicker_dimnames <- function(x)
{
    seed_dimnames <- dimnames(x@seed)
    seed_dimnames[x@picked_dims]  # will return NULL if 'seed_dimnames' is NULL
}

setMethod("dimnames", "SeedDimPicker", .get_SeedDimPicker_dimnames)

.extract_array_from_SeedDimPicker <- function(x, index)
{
    seed_dim <- dim(x@seed)
    seed_index <- rep.int(list(1L), length(seed_dim))
    seed_index[x@picked_dims] <- index
    subseed <- extract_array(x@seed, seed_index)
    dim(subseed) <- dim(subseed)[sort(x@picked_dims)]
    aperm(subseed, perm=rank(x@picked_dims))
}

setMethod("extract_array", "SeedDimPicker", .extract_array_from_SeedDimPicker)

