### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Manipulating an Nindex
###
### An Nindex is a "multidimensional subsetting index" is a list with one
### subscript per dimension in the object to subset. Missing subscripts are
### represented by NULLs. Before it can be used to actually subset an
### array-like object, the NULLs must be replaced with elements of class "name".

### Used in HDF5Array!
expand_Nindex_RangeNSBS <- function(Nindex)
{
    stopifnot(is.list(Nindex))
    RangeNSBS_idx <- which(vapply(Nindex, is, logical(1), "RangeNSBS"))
    Nindex[RangeNSBS_idx] <- lapply(Nindex[RangeNSBS_idx], as.integer)
    Nindex
}

.Nindex2subscripts <- function(Nindex, x)
{
    stopifnot(is.list(Nindex), length(Nindex) == length(dim(x)))

    if (is.array(x))
        Nindex <- expand_Nindex_RangeNSBS(Nindex)

    ## Replace NULLs with list elements of class "name".
    subscripts <- rep.int(alist(foo=), length(Nindex))
    names(subscripts ) <- names(Nindex)
    not_missing_idx <- which(!S4Vectors:::sapply_isNULL(Nindex))
    subscripts[not_missing_idx] <- Nindex[not_missing_idx]

    subscripts
}

subset_by_Nindex <- function(x, Nindex, drop=FALSE)
{
    Nindex <- .Nindex2subscripts(Nindex, x)
    do.call(`[`, c(list(x), Nindex, list(drop=drop)))
}

replace_by_Nindex <- function(x, Nindex, value)
{
    Nindex <- .Nindex2subscripts(Nindex, x)
    do.call(`[<-`, c(list(x), Nindex, list(value=value)))
}

Nindex_as_string <- function(Nindex, dimnames=NULL)
{
    stopifnot(is.list(Nindex))
    missing_idx <- which(S4Vectors:::sapply_isNULL(Nindex))
    Nindex[missing_idx] <- ""
    s <- as.character(Nindex)
    RangeNSBS_idx <- which(vapply(Nindex, is, logical(1), "RangeNSBS"))
    s[RangeNSBS_idx] <- lapply(Nindex[RangeNSBS_idx],
        function(i) paste0(i@subscript, collapse=":")
    )
    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(Nindex) == length(dimnames))
        usename_idx <- which(nzchar(s) &
                             lengths(Nindex) == 1L &
                             lengths(dimnames) != 0L)
        s[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                      Nindex[usename_idx],
                                      SIMPLIFY=FALSE)
    }
    paste0(s, collapse=", ")
}

### Used in HDF5Array!
### Return the lengths of the subscripts in 'Nindex'. The length of a
### missing subscript is the length it would have after expansion.
get_Nindex_lengths <- function(Nindex, dim)
{
    stopifnot(is.list(Nindex), length(Nindex) == length(dim))
    ans <- lengths(Nindex)
    missing_idx <- which(S4Vectors:::sapply_isNULL(Nindex))
    ans[missing_idx] <- dim[missing_idx]
    ans
}

### 'dimnames' must be NULL or a list of the same length as 'Nindex'.
### 'along' must be an integer >= 1 and <= length(Nindex).
get_Nindex_names_along <- function(Nindex, dimnames, along)
{
    stopifnot(is.list(Nindex))
    i <- Nindex[[along]]
    if (is.null(i))
        return(dimnames[[along]])
    names(i)
}

### Used in HDF5Array!
### Return an Nindex ("multidimensional subsetting index").
make_Nindex_from_block_ranges <- function(block_ranges, dim,
                                          expand.RangeNSBS=FALSE)
{
    stopifnot(is(block_ranges, "Ranges"),
              is.integer(dim))
    ndim <- length(dim)
    block_offsets <- start(block_ranges)
    block_dim <- width(block_ranges)
    stopifnot(length(block_ranges) == ndim,
              all(block_offsets >= 1L),
              all(block_offsets + block_dim - 1L <= dim))
    Nindex <- vector(mode="list", length=ndim)
    is_not_missing <- block_dim < dim
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        is_width1 <- block_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
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
    Nindex[expand_idx] <- as.list(block_ranges[expand_idx])
    Nindex
}

### Convert 'Nindex' to a "linear index".
to_linear_index <- function(Nindex, dim)
{
    stopifnot(is.list(Nindex), is.integer(dim), length(Nindex) == length(dim))
    i <- p <- 1L
    for (n in seq_along(Nindex)) {
        idx <- Nindex[[n]]
        d <- dim[[n]]
        if (is.null(idx))
            idx <- seq_len(d)
        i <- rep((idx - 1L) * p, each=length(i)) + i
        p <- p * d
    }
    i
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
    part_idx <- findInterval(idx, breakpoints + 1L) + 1L
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

