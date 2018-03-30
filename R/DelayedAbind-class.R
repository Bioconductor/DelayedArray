### =========================================================================
### DelayedAbind objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("DelayedAbind",
    representation(
        seeds="list",    # List of array-like objects to bind. Each object
                         # is expected to satisfy the "seed contract" i.e.
                         # to support dim(), dimnames(), and extract_array().

        along="integer"  # Single integer indicating the dimension along
                         # which to bind the seeds.
    ),
    prototype(
        seeds=list(new("array")),
        along=1L
    )
)

.validate_DelayedAbind <- function(x)
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

setValidity2("DelayedAbind", .validate_DelayedAbind)

new_DelayedAbind <- function(seeds, along)
{
    seeds <- lapply(seeds, remove_pristine_DelayedArray_wrapping)
    new2("DelayedAbind", seeds=seeds, along=along)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array()
###

.get_DelayedAbind_dim <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dims_along(dims, x@along)
}

setMethod("dim", "DelayedAbind", .get_DelayedAbind_dim)

.get_DelayedAbind_dimnames <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dimnames_along(x@seeds, dims, x@along)
}

setMethod("dimnames", "DelayedAbind", .get_DelayedAbind_dimnames)

.extract_array_from_DelayedAbind <- function(x, index)
{
    i <- index[[x@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        tmp <- lapply(x@seeds, extract_array, index)
        ## Bind the ordinary arrays in 'tmp'.
        ans <- do.call(simple_abind, c(tmp, list(along=x@along)))
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- get_part_index(i, breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[x@along]] <- split_part_idx[[s]]
        extract_array(x@seeds[[s]], index)
    }
    tmp <- lapply(seq_along(x@seeds), FUN)

    ## Bind the ordinary arrays in 'tmp'.
    ans <- do.call(simple_abind, c(tmp, list(along=x@along)))

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector(mode="list", length=length(index))
    Nindex[[x@along]] <- get_rev_index(part_idx)
    subset_by_Nindex(ans, Nindex)
}

setMethod("extract_array", "DelayedAbind", .extract_array_from_DelayedAbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###
### DelayedAbind objects don't support the path() getter/setter for now.
###

setMethod("path", "DelayedAbind",
    function(object, ...)
        stop(wmsg("path() is not supported on a DelayedArray ",
                  "object with multiple leaf seeds at the moment"))
)

setReplaceMethod("path", "DelayedAbind",
    function(object, ..., value)
        stop(wmsg("the path() setter is not supported on a DelayedArray ",
                  "object with multiple leaf seeds at the moment"))
)

