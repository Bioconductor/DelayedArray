### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Using a "multidimensional subscript"
###

### Used in HDF5Array!
### 'subscripts' must be a "multidimensional subscript" i.e. a list with one
### subscript per dimension in 'x'. Missing subscripts are represented by
### list elements of class "name".
subset_by_subscripts <- function(x, subscripts, drop=FALSE)
{
    stopifnot(is.list(subscripts), length(subscripts) == length(dim(x)))
    do.call(`[`, c(list(x), subscripts, list(drop=drop)))
}

expand_missing_subscripts <- function(subscripts, dim)
{
    stopifnot(is.list(subscripts), length(subscripts) == length(dim))
    missing_idx <- which(vapply(subscripts, is.name, logical(1)))
    subscripts[missing_idx] <- lapply(dim[missing_idx], seq_len)
    subscripts
}

### Used in HDF5Array!
### Return the lengths of the subscripts in 'subscripts'. The length of a
### missing subscript is the length it would have after expansion.
### In other words, 'get_subscripts_lengths(subscripts, dim)' is equivalent
### to 'lengths(expand_missing_subscripts(subscripts, dim))' but is more
### efficient because it doesn't actually expand missing subscripts.
get_subscripts_lengths <- function(subscripts, dim)
{
    stopifnot(is.list(subscripts), length(subscripts) == length(dim))
    missing_idx <- which(vapply(subscripts, is.name, logical(1)))
    ans <- lengths(subscripts)
    ans[missing_idx] <- dim[missing_idx]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Translate an index into the whole to an index into the parts
###
### This is .rowidx2rowkeys() from BSgenome/R/OnDiskLongTable-class.R, copied
### here and renamed get_part_index() !
### TODO: Put it somewhere else where it can be shared.
###

.breakpoints2offsets <- function(breakpoints)
{
    breakpoints_len <- length(breakpoints)
    if (breakpoints_len == 0L)
        return(integer(0))
    c(0L, breakpoints[-breakpoints_len])
}

### breakpoints: integer vector of break points that define the parts.
### idx:         index into the whole as an integer vector.
### Return a list of 2 integer vectors parallel to 'idx'. The 1st vector
### contains part numbers and the 2nd vector indices into the parts.
### In addition, if 'breakpoints' has names (part names) then they are
### propagated to the 1st vector.
get_part_index <- function(idx, breakpoints)
{
    part_idx <- findInterval(idx - 1L, breakpoints) + 1L
    names(part_idx) <- names(breakpoints)[part_idx]
    rel_idx <- idx - .breakpoints2offsets(unname(breakpoints))[part_idx]
    list(part_idx, rel_idx)
}

split_part_index <- function(part_index, npart)
{
    ans <- rep.int(list(integer(0)), npart)
    tmp <- split(unname(part_index[[2L]]), part_index[[1L]])
    ans[as.integer(names(tmp))] <- tmp
    ans
}

get_rev_index <- function(part_index)
{
    f <- part_index[[1L]]
    idx <- split(seq_along(f), f)
    idx <- unlist(idx, use.names=FALSE)
    rev_idx <- integer(length(idx))
    rev_idx[idx] <- seq_along(idx)
    rev_idx
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### require_HDF5Array()
###

require_HDF5Array <- function()
{
    if (!requireNamespace("HDF5Array", quietly=TRUE))
        stop("This action requires the HDF5Array package. Please ",
             "install it with:\n\n    library(BiocInstaller)\n",
             "    biocLite(\"HDF5Array\")\n\n  and try again.")
}

