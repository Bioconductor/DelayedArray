### =========================================================================
### DelayedSubset objects
### -------------------------------------------------------------------------
###
### Representation of a delayed multi-dimensional single bracket subsetting
### operation.
###

setClass("DelayedSubset",
    contains="DelayedUnaryOp",
    representation(
        index="list"  # List of subscripts as positive integer vectors,
                      # one per dimension in the input. **Missing** list
                      # elements are allowed and represented by NULLs.
    ),
    prototype(
        index=list(NULL)
    )
)

.validate_DelayedSubset <- function(x)
{
    ## 'index' slot.
    if (length(x@index) != length(dim(x@seed)))
        return("'x@index' must have one list element per dimension in 'x@seed'")
    if (!is.null(names(x@index)))
        return("'x@index' should not have names")
    ok <- lapply(x@index,
              function(i) {is.null(i) || is.integer(i) && is.null(names(i))})
    if (!all(unlist(ok)))
        return(paste0("each list element in 'x@index' must be NULL ",
                      "or an integer vector with no names on it"))
    TRUE
}

setValidity2("DelayedSubset", .validate_DelayedSubset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### 'Nindex' must be a "multidimensional subsetting Nindex" (see
### R/Nindex-utils.R in the S4Arrays package) or NULL.
new_DelayedSubset <- function(seed=new("array"), Nindex=NULL)
{
    index <- S4Arrays:::normalize_Nindex(Nindex, seed)
    new2("DelayedSubset", seed=seed, index=index)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_noop() method
###

setMethod("is_noop", "DelayedSubset",
    function(x) all(S4Vectors:::sapply_isNULL(x@index))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedSubset

.DelayedSubset_summary <- function(object) "Subset"

summary.DelayedSubset <-
    function(object, ...) .DelayedSubset_summary(object, ...)

setMethod("summary", "DelayedSubset", summary.DelayedSubset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###

setMethod("dim", "DelayedSubset",
    function(x) S4Arrays:::get_Nindex_lengths(x@index, dim(x@seed))
)

setMethod("dimnames", "DelayedSubset",
    function(x) S4Arrays:::subset_dimnames_by_Nindex(dimnames(x@seed), x@index)
)

subset_DelayedSubset <- function(x, index)
{
    stopifnot(is(x, "DelayedSubset"))
    x_ndim <- length(x@index)
    stopifnot(is.list(index), length(index) == x_ndim)
    seed_dim <- dim(x@seed)
    ## Would mapply() be faster here?
    x@index <- lapply(seq_len(x_ndim),
        function(along) {
            i0 <- x@index[[along]]
            i <- index[[along]]
            if (is.null(i))
                return(i0)
            if (is.null(i0))
                return(i)
            ans <- i0[i]
            if (isSequence(ans, of.length=seed_dim[[along]]))
                return(NULL)
            ans
        })
    x
}

setMethod("extract_array", "DelayedSubset",
    function(x, index)
    {
        x2 <- subset_DelayedSubset(x, index)
        extract_array(x2@seed, x2@index)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedSubset",
    function(x)
    {
        if (!is_sparse(x@seed))
            return(FALSE)
        ## Duplicates in x@index break structural sparsity.
        !any(vapply(x@index, anyDuplicated,
                    integer(1), USE.NAMES=FALSE))
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedSubset",
    function(x, index)
    {
        x2 <- subset_DelayedSubset(x, index)
        ## Assuming that the caller respected "extract_sparse_array() Terms
        ## of Use" (see SparseArraySeed-class.R), 'is_sparse(x)' should be
        ## TRUE and the subscripts in 'index' should not contain duplicates.
        ## This in turn means that the subscripts in 'x2@index' should not
        ## contain duplicates either so the call below should also respect
        ## "extract_sparse_array() Terms of Use".
        extract_sparse_array(x2@seed, x2@index)
    }
)

