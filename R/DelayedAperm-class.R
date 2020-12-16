### =========================================================================
### DelayedAperm objects
### -------------------------------------------------------------------------
###
### Representation of a delayed "extended aperm()" operation, that is, a
### delayed aperm() that can drop and/or add **ineffective** dimensions.
### Note that since only **ineffective** dimensions (i.e. dimensions with
### an extent of 1) can be dropped or added, the length of the output array
### is guaranteed to be the same as the length of the input array.
###

setClass("DelayedAperm",
    contains="DelayedUnaryOp",
    representation(
        perm="integer"  # Index into 'dim(seed)' describing the
                        # **rearrangement** of the dimensions i.e. which
                        # dimensions of the input to keep and in which order.
                        # Only ineffective dimensions can be dropped. Note
                        # that NAs are allowed and indicate the addition of
                        # an ineffective dimension. For example if 'dim(seed)'
                        # is c(20, 1, 15, 2, 1) then a DelayedAperm object
                        # where 'perm' is set to c(NA, NA, 3, 1, NA, 4, 5)
                        # represents an operation that returns an array with
                        # dimensions c(1, 1, 15, 20, 1, 2, 1).
    ),
    prototype(
        perm=1L
    )
)

.validate_DelayedAperm <- function(x)
{
    ## 'perm' slot.
    msg <- validate_perm(x@perm, dim(x@seed))
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("DelayedAperm", .validate_DelayedAperm)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

new_DelayedAperm <- function(seed=new("array"), perm=NULL)
{
    perm <- normarg_perm(perm, dim(seed))
    new2("DelayedAperm", seed=seed, perm=perm)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_noop() method
###

setMethod("is_noop", "DelayedAperm",
    function(x) isSequence(x@perm, length(dim(x@seed)))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedAperm

.DelayedAperm_summary <- function(object)
{
    perm <- as.character(object@perm)
    if (length(perm) >= 2L)
        perm <- sprintf("c(%s)", paste0(perm, collapse=","))
    sprintf("Aperm (perm=%s)", perm)
}

summary.DelayedAperm <-
    function(object, ...) .DelayedAperm_summary(object, ...)

setMethod("summary", "DelayedAperm", summary.DelayedAperm)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###

.get_DelayedAperm_dim <- function(x)
{
    seed_dim <- dim(x@seed)
    ans <- seed_dim[x@perm]
    ans[is.na(x@perm)] <- 1L
    ans
}

setMethod("dim", "DelayedAperm", .get_DelayedAperm_dim)

.get_DelayedAperm_dimnames <- function(x)
{
    seed_dimnames <- dimnames(x@seed)
    if (is.null(seed_dimnames))
        return(NULL)
    simplify_NULL_dimnames(seed_dimnames[x@perm])
}

setMethod("dimnames", "DelayedAperm", .get_DelayedAperm_dimnames)

project_index_on_seed <- function(index, x)
{
    stopifnot(is(x, "DelayedAperm"),
              is.list(index),
              length(index) == length(x@perm))
    nonNA_idx <- which(!is.na(x@perm))
    perm0 <- x@perm[nonNA_idx]
    index0 <- index[nonNA_idx]
    seed_dim <- dim(x@seed)
    seed_index <- vector("list", length=length(seed_dim))
    seed_index[perm0] <- index0
    seed_index
}

.extract_array_from_DelayedAperm <- function(x, index)
{
    seed_index <- project_index_on_seed(index, x)
    a <- extract_array(x@seed, seed_index)
    a <- aperm2(a, x@perm)
    index[!is.na(x@perm)] <- list(NULL)
    subset_by_Nindex(a, index)
}

setMethod("extract_array", "DelayedAperm",
    .extract_array_from_DelayedAperm
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedAperm", function(x) is_sparse(x@seed))

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedAperm",
    function(x, index)
    {
        seed_index <- project_index_on_seed(index, x)
        sas <- extract_sparse_array(x@seed, seed_index)
        sas <- aperm(sas, x@perm)
        index[!is.na(x@perm)] <- list(NULL)
        extract_sparse_array(sas, index)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Backward compatibility with DelayedArray < 0.5.24
###
### In DelayedArray 0.5.24 the SeedDimPicker class got renamed DelayedAperm.
### DelayedArray objects serialized with DelayedArray < 0.5.24 might contain
### SeedDimPicker instances nested in their "seed" slot so we need to keep
### the class around for now.
###

setClass("SeedDimPicker", contains="DelayedAperm")

setMethod("updateObject", "SeedDimPicker",
    function(object, ..., verbose=FALSE)
    {
        object <- new2("DelayedAperm", seed=object@seed,
                                       perm=object@dim_combination)
        callNextMethod()
    }
)

