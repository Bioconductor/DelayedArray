### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### TODO: Move this to S4Vectors and implment it in C. Does not need to
### create 'seq_len(d)' and will be able to do early bailout.
is_sequence <- function(x, length)
{
    stopifnot(is.integer(x), isSingleInteger(length))
    length(x) == length && identical(x, seq_len(length))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 wrappers to dim<- and dimnames<- that try to avoid unnecessary copies
### of 'x'
###

set_dim <- function(x, value)
{
    if (!identical(dim(x), value))
        dim(x) <- value
    x
}

set_dimnames <- function(x, value)
{
    if (!identical(dimnames(x), value))
        dimnames(x) <- value
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_dimnames()
###

simplify_NULL_dimnames <- function(dimnames)
{
    if (all(S4Vectors:::sapply_isNULL(dimnames)))
        return(NULL)
    dimnames
}

### Like for extract_array(), 'index' is expected to be an unnamed list of
### subscripts as positive integer vectors, one vector per dimension in 'x'.
### *Missing* list elements are allowed and represented by NULLs.
### See extract_array.R for more information.
subset_dimnames <- function(dimnames, index)
{
    stopifnot(is.list(index))
    if (is.null(dimnames))
        return(NULL)
    ndim <- length(index)
    stopifnot(is.list(dimnames), length(dimnames) == ndim)
    ## Would mapply() be faster here?
    ans <- lapply(seq_len(ndim),
                  function(along) {
                      dn <- dimnames[[along]]
                      i <- index[[along]]
                      if (is.null(dn) || is.null(i))
                          return(dn)
                      dn[i]
                  })
    simplify_NULL_dimnames(ans)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Manipulating an Nindex
###
### An Nindex is a "multidimensional subsetting index". It's represented as a
### list with one subscript per dimension in the array-like object to subset.
### NULL list elements in it are interpreted as missing subscripts, that is, as
### subscripts that run along the full extend of the corresponding dimension.
### Before an Nindex can be used in a call to `[`, `[<-`, `[[` or `[[<-`, the
### NULL list elements must be replaced with object of class "name".
###

### For use in "[", "[<-", "[[", or "[[<-" methods to extract the user
### supplied subscripts as an Nindex. NULL subscripts are replace with
### integer(0). Missing subscripts are set to NULL.
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
        Nindex <- Nindex[!(argnames %in% c("drop", "exact", "value"))]
    if (length(Nindex) == 1L && is.null(Nindex[[1L]]))
        Nindex <- Nindex[0L]
    Nindex
}

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
    subscripts <- rep.int(list(quote(expr=)), length(Nindex))
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

### Return the modified array.
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_array_objects()
###

### 'objects' must be a list of array-like objects that support as.vector().
combine_array_objects <- function(objects)
{
    if (!is.list(objects))
        stop("'objects' must be a list")
    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    if (length(objects) == 0L)
        return(NULL)
    unlist(lapply(objects, as.vector), recursive=FALSE, use.names=FALSE)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used in validity methods
###

### A modified version of S4Vectors::wmsg() that is better suited for use
### by validity methods.
### TODO: Put this in S4Vectors next to wmsg(). Would probably need a better
### name.
wmsg2 <- function(...)
    paste0("\n    ",
           paste0(strwrap(paste0(c(...), collapse="")), collapse="\n    "))

validate_dim_slot <- function(x, slotname="dim")
{
    x_dim <- slot(x, slotname)
    if (!is.integer(x_dim))
        return(wmsg2(sprintf("'%s' slot must be an integer vector", slotname)))
    if (length(x_dim) == 0L)
        return(wmsg2(sprintf("'%s' slot cannot be empty", slotname)))
    if (S4Vectors:::anyMissingOrOutside(x_dim, 0L))
        return(wmsg2(sprintf("'%s' slot cannot contain negative or NA values",
                             slotname)))
    TRUE
}

validate_dimnames_slot <- function(x, dim, slotname="dimnames")
{
    x_dimnames <- slot(x, slotname)
    if (!is.list(x_dimnames))
        return(wmsg2(sprintf("'%s' slot must be a list", slotname)))
    if (length(x_dimnames) != length(dim))
        return(wmsg2(sprintf("'%s' slot must have ", slotname),
                     "one list element per dimension in the object"))
    ok <- vapply(seq_along(dim),
                 function(along) {
                   dn <- x_dimnames[[along]]
                   if (is.null(dn))
                       return(TRUE)
                   is.character(dn) && length(dn) == dim[[along]]
                 },
                 logical(1),
                 USE.NAMES=FALSE)
    if (!all(ok))
        return(wmsg2(sprintf("each list element in '%s' slot ", slotname),
                     "must be NULL or a character vector along ",
                     "the corresponding dimension in the object"))
    TRUE
}

