### =========================================================================
### DelayedSubset objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("DelayedSubset",
    representation(
        seed="ANY",   # An array-like object expected to satisfy the
                      # "seed contract" i.e. to support dim(), dimnames(),
                      # and extract_array().

        index="list"  # List (possibly named) of subscripts as
                      # positive integer vectors, one vector per
                      # seed dimension. *Missing* list elements
                      # are allowed and represented by NULLs.
    ),
    prototype(
        seed=new("array"),
        index=list(NULL)
    )
)

.validate_DelayedSubset <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)
    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))
    ## 'index' slot.
    if (length(x@index) != seed_ndim)
        return(wmsg2("'x@index' must have one list element per dimension ",
                     "in 'x@seed'"))
    if (!all(S4Vectors:::sapply_isNULL(x@index) |
             vapply(x@index, is.integer, logical(1), USE.NAMES=FALSE)))
        return(wmsg2("every list element in 'x@index' must be either NULL ",
                     "or an integer vector"))
    TRUE
}

setValidity2("DelayedSubset", .validate_DelayedSubset)

new_DelayedSubset <- function(seed=new("array"), index=list(NULL))
{
    new2("DelayedSubset", seed=seed, index=index)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array()
###

.get_DelayedSubset_dim <- function(x) get_Nindex_lengths(x@index, dim(x@seed))

setMethod("dim", "DelayedSubset", .get_DelayedSubset_dim)

.get_DelayedSubset_dimnames <- function(x)
{
    x_seed_dimnames <- dimnames(x@seed)
    ans <- lapply(seq_along(x@index),
                  get_Nindex_names_along,
                      Nindex=x@index,
                      dimnames=x_seed_dimnames)
    if (all(S4Vectors:::sapply_isNULL(ans)))
        return(NULL)
    ans
}

setMethod("dimnames", "DelayedSubset", .get_DelayedSubset_dimnames)

.extract_array_from_DelayedSubset <- function(x, index)
{
    x_seed_dim <- dim(x@seed)
    stopifnot(is.list(index), length(index) == length(x_seed_dim))
    index2 <- lapply(seq_along(x@index),
                     function(along) {
                         i1 <- x@index[[along]]
                         i2 <- index[[along]]
                         if (is.null(i2))
                             return(i1)
                         if (is.null(i1))
                             i1 <- seq_len(x_seed_dim[[along]])  # expand 'i1'
                         i1[i2]
                     })
    extract_array(x@seed, index2)
}

setMethod("extract_array", "DelayedSubset", .extract_array_from_DelayedSubset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seed() getter/setter
###

setMethod("seed", "DelayedSubset", function(x) x@seed)

setReplaceMethod("seed", "DelayedSubset",
    function(x, value)
    {
        x@seed <- normalize_seed_replacement_value(value, seed(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### The path of a DelayedSubset object is the path of its seed so path()
### will work only on a DelayedSubset object with a seed that supports
### path().
### For example it will work if the seed is an on-disk object (e.g. an
### HDF5ArraySeed object) but not if it's an in-memory object (e.g. an
### ordinary array or RleArraySeed object).
setMethod("path", "DelayedSubset",
    function(object, ...) path(seed(object), ...)
)

### The path() setter sets the supplied path on the seed of the
### DelayedSubset object so it will work out-of-the-box on any
### DelayedSubset object with a seed that supports the path() setter.
### For example it will work if the seed is an HDF5ArraySeed object.
setReplaceMethod("path", "DelayedSubset",
    function(object, ..., value)
    {
        path(seed(object), ...) <- value
        object
    }
)

