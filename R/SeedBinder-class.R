### =========================================================================
### SeedBinder objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("SeedBinder",
    representation(
        seeds="list",    # List of array-like objects to bind. Each object
                         # is expected to satisfy the "seed contract" i.e.
                         # to support dim(), dimnames(), and
                         # subset_seed_as_array().

        along="integer"  # Single integer indicating the dimension along
                         # which to bind the seeds.
    ),
    prototype(
        seeds=list(new("array")),
        along=1L
    )
)

.validate_SeedBinder <- function(x)
{
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))
    if (!(isSingleInteger(x@along) && x@along > 0L))
        return(wmsg2("'x@along' must be a single positive integer"))
    dims <- get_dims_to_bind(x@seeds, x@along)
    if (is.character(dims))
        return(wmsg2(dims))
    TRUE
}

setValidity2("SeedBinder", .validate_SeedBinder)

new_SeedBinder <- function(seeds, along)
{
    seeds <- lapply(seeds, remove_pristine_DelayedArray_wrapping)
    new2("SeedBinder", seeds=seeds, along=along)
}

### Implement the "seed contract" i.e. dim(), dimnames(), and
### subset_seed_as_array().

.get_SeedBinder_dim <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dims_along(dims, x@along)
}

setMethod("dim", "SeedBinder", .get_SeedBinder_dim)

.get_SeedBinder_dimnames <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dimnames_along(x@seeds, dims, x@along)
}

setMethod("dimnames", "SeedBinder", .get_SeedBinder_dimnames)

.subset_SeedBinder_as_array <- function(seed, index)
{
    i <- index[[seed@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        tmp <- lapply(seed@seeds, subset_seed_as_array, index)
        ## Bind the ordinary arrays in 'tmp'.
        ans <- do.call(simple_abind, c(tmp, list(along=seed@along)))
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- get_dims_to_bind(seed@seeds, seed@along)
    breakpoints <- cumsum(dims[seed@along, ])
    part_idx <- get_part_index(i, breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[seed@along]] <- split_part_idx[[s]]
        subset_seed_as_array(seed@seeds[[s]], index)
    }
    tmp <- lapply(seq_along(seed@seeds), FUN)

    ## Bind the ordinary arrays in 'tmp'.
    ans <- do.call(simple_abind, c(tmp, list(along=seed@along)))

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector(mode="list", length=length(index))
    Nindex[[seed@along]] <- get_rev_index(part_idx)
    subset_by_Nindex(ans, Nindex)
}

setMethod("subset_seed_as_array", "SeedBinder", .subset_SeedBinder_as_array)

