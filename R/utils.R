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
expand_RangeNSBS_subscripts <- function(subscripts)
{
    stopifnot(is.list(subscripts))
    RangeNSBS_idx <- which(vapply(subscripts, is, logical(1), "RangeNSBS"))
    subscripts[RangeNSBS_idx] <- lapply(subscripts[RangeNSBS_idx], as.integer)
    subscripts
}

### Used in HDF5Array!
### Return a "multidimensional subscript" i.e. a list with one subscript per
### dimension in the original array. Missing subscripts are represented by
### list elements of class "name".
make_subscripts_from_ranges <- function(block_ranges, dim,
                                        expand.RangeNSBS=FALSE)
{
    stopifnot(is(block_ranges, "Ranges"),
              is.integer(dim),
              length(block_ranges) == length(dim))
    block_offsets <- start(block_ranges)
    block_dim <- width(block_ranges)
    stopifnot(all(block_dim <= dim))
    ndim <- length(dim)
    subscripts <- rep.int(alist(foo=), ndim)
    is_not_missing <- block_dim != dim
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        is_width1 <- block_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        subscripts[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
            function(i) {
                start <- block_offsets[[i]]
                end <- start + block_dim[[i]] - 1L
                upper_bound <- dim[[i]]
                new2("RangeNSBS", subscript=c(start, end),
                                  upper_bound=upper_bound,
                                  check=FALSE)
            }
        )
    }
    subscripts[expand_idx] <- as.list(block_ranges[expand_idx])
    subscripts
}

### Used in HDF5Array!
### 'subscripts' must be a "multidimensional subscript" i.e. a list with one
### subscript per dimension in 'x'. Missing subscripts are represented by
### list elements of class "name".
subset_by_subscripts <- function(x, subscripts, drop=FALSE)
{
    stopifnot(is.list(subscripts), length(subscripts) == length(dim(x)))
    if (is.array(x))
        subscripts <- expand_RangeNSBS_subscripts(subscripts)
    do.call(`[`, c(list(x), subscripts, list(drop=drop)))
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

subscripts_as_string <- function(subscripts, dimnames=NULL)
{
    stopifnot(is.list(subscripts))
    s <- as.character(subscripts)
    RangeNSBS_idx <- which(vapply(subscripts, is, logical(1), "RangeNSBS"))
    s[RangeNSBS_idx] <- lapply(subscripts[RangeNSBS_idx],
        function(i) paste0(i@subscript, collapse=":")
    )
    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(subscripts) == length(dimnames))
        usename_idx <- which(nzchar(s) &
                             lengths(subscripts) == 1L &
                             lengths(dimnames) != 0L)
        s[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                      subscripts[usename_idx],
                                      SIMPLIFY=FALSE)
    }
    paste0(s, collapse=", ")
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

