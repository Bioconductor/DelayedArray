### =========================================================================
### Manipulating an Nindex
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###
### An Nindex is a "multidimensional subsetting index". It's represented as a
### list with one subscript per dimension in the array-like object to subset.
### NULL list elements in it are interpreted as missing subscripts, that is, as
### subscripts that run along the full extend of the corresponding dimension.
### Before an Nindex can be used in a call to `[`, `[<-`, `[[` or `[[<-`, the
### NULL list elements must be replaced with objects of class "name".
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Normalization of an Nindex
###

normalizeSingleBracketSubscript2 <- function(i, x_len, x_names=NULL)
{
    ## We support subsetting by an array-like subscript but only if the
    ## subscript is mono-dimensional, in which case we call as.vector() on
    ## it. This will possibly trigger its realization e.g. if it's a
    ## DelayedArray object.
    i_dim <- dim(i)
    if (!is.null(i_dim)) {
        if (length(i_dim) != 1L)
            stop(wmsg("subsetting a DelayedArray object with an array-like ",
                      "subscript is only supported if the subscript has a ",
                      "single dimension"))
        i <- as.vector(i)
    }
    ## We create an artificial object 'x' of length 'x_len' with 'x_names' on
    ## it. normalizeSingleBracketSubscript() will only look at its length and
    ## names so what the object really is doesn't matter. Hence we make it
    ## with the smallest possible memory footprint.
    ## TODO: Change the signature of normalizeSingleBracketSubscript() in
    ## S4Vectors to take 'x_len' and 'x_names' instead of 'x' so we won't
    ## have to use this kind of trick.
    if (is.null(x_names)) {
        x <- Rle(0L, x_len)
    } else {
        x <- setNames(raw(x_len), x_names)
    }
    normalizeSingleBracketSubscript(i, x)
}

### Normalize 'Nindex' i.e. check and turn each non-NULL list element
### into a positive integer vector that is a valid subscript along the
### corresponding dimension in 'x'.
normalizeNindex <- function(Nindex, x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    x_ndim <- length(x_dim)
    if (is.null(Nindex))
        return(vector("list", length=x_ndim))
    if (!(is.list(Nindex) && length(Nindex) == x_ndim))
        stop(wmsg("'Nindex' must be a list with one ",
                  "list element per dimension in 'x'"))
    x_dimnames <- dimnames(x)
    lapply(seq_len(x_ndim),
        function(along) {
            subscript <- Nindex[[along]]
            if (is.null(subscript))
                return(NULL)
            d <- x_dim[[along]]
            i <- normalizeSingleBracketSubscript2(subscript, d,
                                                  x_dimnames[[along]])
            if (isSequence(i, of.length=d))
                return(NULL)
            i
        })
}

### Assume 'Nindex' is normalized (see above).
### Return a logical vector with one logical per dimension indicating
### whether the corresponding subscript in 'Nindex' reaches all valid
### positions along the dimension.
subscript_has_nogap <- function(Nindex, dim)
{
    stopifnot(is.list(Nindex), length(Nindex) == length(dim))
    vapply(seq_along(Nindex),
        function(along) {
            Li <- Nindex[[along]]
            if (is.null(Li))
                return(TRUE)
            d <- dim[[along]]
            if (length(Li) < d)
                return(FALSE)
            hits <- logical(d)
            hits[Li] <- TRUE
            all(hits)
        },
        logical(1L),
        USE.NAMES=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other Nindex utilities
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
    expand_idx <- which(vapply(Nindex, is, logical(1), "RangeNSBS"))
    if (length(expand_idx) != 0L)
        Nindex[expand_idx] <- lapply(Nindex[expand_idx], as.integer)
    Nindex
}

.make_subscripts_from_Nindex <- function(Nindex, x)
{
    stopifnot(is.list(Nindex), length(Nindex) == length(dim(x)))

    if (is.array(x))
        Nindex <- expand_Nindex_RangeNSBS(Nindex)

    ## Replace NULLs with list elements of class "name".
    subscripts <- rep.int(list(quote(expr=)), length(Nindex))
    names(subscripts) <- names(Nindex)
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

subset_dimnames_by_Nindex <- function(dimnames, Nindex)
{
    stopifnot(is.list(Nindex))
    if (is.null(dimnames))
        return(NULL)
    ndim <- length(Nindex)
    stopifnot(is.list(dimnames), length(dimnames) == ndim)
    ## Would mapply() be faster here?
    ans <- lapply(setNames(seq_len(ndim), names(dimnames)),
                  function(along) {
                      dn <- dimnames[[along]]
                      i <- Nindex[[along]]
                      if (is.null(dn) || is.null(i))
                          return(dn)
                      extractROWS(dn, i)
                  })
    simplify_NULL_dimnames(ans)
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

