### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A modified version of S4Vectors::wmsg() that is better suited for use
### by validity methods.
### TODO: Put this in S4Vectors next to wmsg(). Would probably need a better
### name.

wmsg2 <- function(...)
    paste0("\n    ",
           paste0(strwrap(paste0(c(...), collapse="")), collapse="\n    "))


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

.make_subscripts_from_Nindex <- function(Nindex, x)
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
    subscripts <- .make_subscripts_from_Nindex(Nindex, x)
    do.call(`[`, c(list(x), subscripts, list(drop=drop)))
}

replace_by_Nindex <- function(x, Nindex, value)
{
    subscripts <- .make_subscripts_from_Nindex(Nindex, x)
    do.call(`[<-`, c(list(x), subscripts, list(value=value)))
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

### Convert 'Nindex' to a "linear index".
### Return the "linear index" as an integer vector if prod(dim) <=
### .Machine$integer.max, otherwise as a vector of doubles.
to_linear_index <- function(Nindex, dim)
{
    stopifnot(is.list(Nindex), is.integer(dim), length(Nindex) == length(dim))
    if (prod(dim) <= .Machine$integer.max) {
        ans <- p <- 1L
    } else {
        ans <- p <- 1
    }
    for (along in seq_along(Nindex)) {
        d <- dim[[along]]
        i <- Nindex[[along]]
        if (is.null(i))
            i <- seq_len(d)
        ans <- rep((i - 1L) * p, each=length(ans)) + ans
        p <- p * d
    }
    ans
}

### For use in "[" or "[[" methods to extract the user supplied subscripts as
### an Nindex. NULL subscripts are replace with integer(0). Missing subscripts
### are set to NULL.
extract_Nindex_from_syscall <- function(call, eframe)
{
    Nindex <- lapply(seq_len(length(call) - 2L),
        function(i) {
            subscript <- call[[2L + i]]
            if (missing(subscript))
                return(NULL)
            subscript <- eval(subscript, envir=eframe, enclos=eframe)
            if (is.null(subscript))
                return(integer(0))
            subscript
        }
    )
    argnames <- tail(names(call), n=-2L)
    if (!is.null(argnames))
        Nindex <- Nindex[!(argnames %in% c("drop", "exact"))]
    if (length(Nindex) == 1L && is.null(Nindex[[1L]]))
        Nindex <- Nindex[0L]
    Nindex
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

