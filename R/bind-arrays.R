### =========================================================================
### Bind arrays with an arbitrary number of dimensions along an arbitrary
### dimension
### -------------------------------------------------------------------------


### Return a matrix with one row per dim and one column per object if the
### objects are bindable. Otherwise return a character vector describing why
### the objects are not bindable. This design allows the function to be used
### in the context of a validity method.
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
        return(c("all the objects to bind must have ",
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
### NOT exported but used in the HDF5Array package.
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
### NOT exported but used in the HDF5Array package.
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
### 'block_lens' is assumed to be an integer vector parallel to 'objects'
### specifying the block length for each object in 'objects'. In addition the
### length of 'object[[i]]' must be 'nblock * block_lens[[i]]' ('nblock' is
### the same for all the objects).
.intertwine_blocks <- function(objects, block_lens)
{
    x0 <- unlist(lapply(objects, `[`, 0L), use.names=FALSE)
    objects_lens <- lengths(objects)
    if (all(objects_lens == 0L))
        return(x0)

    objects <- S4Vectors:::prepare_objects_to_bind(x0, objects)

    nblock <- objects_lens %/% block_lens
    nblock <- unique(nblock[!is.na(nblock)])
    stopifnot(length(nblock) == 1L)  # sanity check

    .Call("abind", objects, nblock, PACKAGE="DelayedArray")
}

### A stripped-down version of abind::abind().
### Some differences:
###   (a) Treatment of dimnames: simple_abind() treatment of dimnames is
###       consistent with base::rbind() and base::cbind(). This is not the
###       case for abind::abind() which does some strange things with the
###       dimnames.
###   (b) Performance: simple_abind() has a little bit more overhead than
###       abind::abind(). This makes it slower on small objects. However it
###       tends to be slightly faster on big objects.
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
    block_lens <- dims[along, ]
    for (n in seq_len(along - 1L))
        block_lens <- block_lens * dims[n, ]
    ans <- .intertwine_blocks(objects, block_lens)

    ## Set the dim.
    ans <- set_dim(ans, combine_dims_along(dims, along))

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

