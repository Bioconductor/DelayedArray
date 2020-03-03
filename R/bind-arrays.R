### =========================================================================
### Bind arrays with an arbitrary number of dimensions along an arbitrary
### dimension
### -------------------------------------------------------------------------


### Return a matrix with one row per dim and one column per object if the
### objects are "bindable". Otherwise return a string describing why they
### are not. This design allows the function to be used in the context of
### a validity method.
get_dims_to_bind <- function(objects, along)
{
    if (!(isSingleInteger(along) && along >= 1L))
        stop("'along' must be a single positive integer")
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    ndim <- ndims[[1L]]
    if (ndim < along)
        stop(wmsg("the array-like objects to bind must have at least ",
                  along, " dimensions for this binding operation"))
    if (!all(ndims == ndim))
        return(paste0("all the objects to bind must have ",
                      "the same number of dimensions"))
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return("the objects to bind have no dimensions")
    dims <- matrix(tmp, nrow=ndim)
    tmp <- dims[-along, , drop=FALSE]
    if (!all(tmp == tmp[ , 1L]))
        return("the objects to bind have incompatible dimensions")
    dims
}

### Combine the dims the rbind/cbind way.
combine_dims_along <- function(dims, along)
{
    stopifnot(is.matrix(dims),
              isSingleInteger(along), along >= 1L, along <= nrow(dims))
    ans_dim <- dims[ , 1L]
    ans_dim[[along]] <- sum(dims[along, ])
    ans_dim
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine the dimnames of a list of array-like objects
###

### Assume all the arrays in 'objects' have the same number of dimensions.
combine_dimnames <- function(objects)
{
    lapply(seq_along(dim(objects[[1L]])),
        function(n) {
            for (x in objects) {
                dn <- dimnames(x)[[n]]
                if (!is.null(dn))
                    return(dn)
            }
            NULL
        })
}

### Combine the dimnames the rbind/cbind way.
combine_dimnames_along <- function(objects, dims, along)
{
    stopifnot(is.matrix(dims),
              isSingleInteger(along), along >= 1L, along <= nrow(dims))
    dimnames <- combine_dimnames(objects)
    along_names <- lapply(objects, function(x) dimnames(x)[[along]])
    along_names_lens <- lengths(along_names)
    if (any(along_names_lens != 0L)) {
        fix_idx <- which(along_names_lens != dims[along, ])
        along_names[fix_idx] <- lapply(dims[along, fix_idx], character)
    }
    along_names <- unlist(along_names, use.names=FALSE)
    if (!is.null(along_names))
        dimnames[[along]] <- along_names
    simplify_NULL_dimnames(dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simple_abind()
###

### 'objects' is assumed to be a list of vector-like objects.
### 'nblock' is assumed to be a single integer value (stored as a numeric)
### that is a common divisor of the object lengths.
.intertwine_blocks <- function(objects, nblock, ans_dim)
{
    x0 <- unlist(lapply(objects, `[`, 0L), recursive=FALSE, use.names=FALSE)
    objects_lens <- lengths(objects)
    if (all(objects_lens == 0L))
        return(set_dim(x0, ans_dim))

    idx <- which(vapply(objects,
        function(object) { typeof(object) != typeof(x0) },
        logical(1),
        USE.NAMES=FALSE))
    if (length(idx) != 0L)
        objects[idx] <- lapply(objects[idx], `storage.mode<-`, typeof(x0))

    .Call2("C_abind", objects, nblock, ans_dim, PACKAGE="DelayedArray")
}

### A stripped-down version of abind::abind().
### Some differences:
###   (a) Treatment of dimnames: simple_abind() treatment of dimnames is
###       consistent with base::rbind() and base::cbind(). This is not the
###       case for abind::abind() which does some strange things with the
###       dimnames.
###   (b) Performance: simple_abind() is much faster than abind::abind()
###       (between 3x and 15x). Also note that in the 'along=1L' and 'along=2L'
###       cases, it's generally as fast (and most of the time faster) than
###       base::rbind() and base::cbind().
###       For example, with 'x <- matrix(1:30000000, nrow=5000)',
###       'simple_abind(m, m, m, along=1L)' is 14x faster than
###       'abind::abind(m, m, m, along=1L)' and 11x faster than
###       'base::rbind(m, m, m)'.
###   (c) abind::abind() is broken on matrices of type "list".
simple_abind <- function(..., along)
{
    objects <- S4Vectors:::delete_NULLs(list(...))
    if (length(objects) == 0L)
        return(NULL)

    ## Check dim compatibility.
    dims <- get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Perform the binding.
    nblock <- prod(dims[-seq_len(along), 1L])  # numeric that can be >
                                               # .Machine$integer.max
    ans <- .intertwine_blocks(objects, nblock, combine_dims_along(dims, along))

    ## Combine and set the dimnames.
    set_dimnames(ans, combine_dimnames_along(objects, dims, along))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bind arrays along their 1st or 2nd dimension
###

setGeneric("arbind", function(...) standardGeneric("arbind"))
setGeneric("acbind", function(...) standardGeneric("acbind"))

setMethod("arbind", "array", function(...) simple_abind(..., along=1L))
setMethod("acbind", "array", function(...) simple_abind(..., along=2L))

