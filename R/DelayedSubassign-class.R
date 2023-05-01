### =========================================================================
### DelayedSubassign objects
### -------------------------------------------------------------------------
###
### Representation of a delayed multi-dimensional single bracket
### subassignment.
###

### Even though strictly speaking DelayedSubassign nodes are binary nodes
### (subassigment operates on 2 array-like objects, the "left value" and the
### "right value"), it turns out to be more convenient (and natural) to treat
### them as unary nodes (e.g. in nseed() and seed()). This is why we make
### DelayedSubassign extend DelayedUnaryOp (via DelayedUnaryIsoOp).
setClass("DelayedSubassign",
    contains="DelayedUnaryIsoOp",
    representation(
        Lindex="list",    # The "left index". List of subscripts as positive
                          # integer vectors, one per dimension in the input.
                          # **Missing** list elements are allowed and
                          # represented by NULLs.
                          # Allowed to contain duplicates BUT NO NAs when the
                          # "Rvalue" slot is an ordinary vector (atomic or
                          # list) of length 1.
                          # Allowed to contain NAs BUT NO DUPLICATES when the
                          # "Rvalue" slot is an array-like object.

        Rvalue="ANY",     # The "right value" i.e. the array-like object on the
                          # right side of the subassignment. Expected to comply
                          # with the "seed contract". Alternatively, it can be
                          # an ordinary vector (atomic or list) of length 1.

        .nogap="logical"  # One logical per dimension in the input indicating
                          # whether the corresponding subscript in the "left
                          # index" reaches all valid positions along the
                          # seed dimension associated with it.
    ),
    prototype(
        Lindex=list(NULL),
        Rvalue=NA,
        .nogap=TRUE
    )
)

.validate_DelayedSubassign <- function(x)
{
    ## TODO!
    TRUE
}

setValidity2("DelayedSubassign", .validate_DelayedSubassign)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_Rvalue <- function(Rvalue, selection_dim)
{
    Rvalue_dim <- dim(Rvalue)
    if (is.null(Rvalue_dim) && !is.vector(Rvalue))
        stop(wmsg("replacement value must be an array-like object ",
                  "or an ordinary vector"))
    ## 'Rvalue' is an array-like object or an ordinary vector (atomic or list).
    if (length(Rvalue) != prod(selection_dim))
        stop(wmsg("length of replacement value must equal the number ",
                  "of array elements to replace"))
    if (is.null(Rvalue_dim)) {
        ## 'x@Rvalue' is an ordinary vector (atomic or list).
        dim(Rvalue) <- selection_dim
        return(Rvalue)
    }
    same_dims <- function(dim1, dim2) length(dim1) == length(dim2) &&
                                      all(dim1 == dim2)
    if (same_dims(Rvalue_dim, selection_dim))
        return(Rvalue)
    ## We're going to reshape 'Rvalue' but only if its effective dimensions
    ## are the same as the effective dimensions of the selection.
    Rvalue_effdim <- Rvalue_dim[Rvalue_dim != 1L]
    selection_effdim <- selection_dim[selection_dim != 1L]
    if (!same_dims(Rvalue_effdim, selection_effdim))
        stop(wmsg("dimensions of replacement value are incompatible ",
                  "with the selection of the subassignment"))
    dim(Rvalue) <- selection_dim
    Rvalue
}

### 'Nindex' must be a "multidimensional subsetting Nindex" (see
### R/Nindex-utils.R in the S4Arrays package) or NULL.
new_DelayedSubassign <- function(seed=new("array"), Nindex=NULL, Rvalue=NA)
{
    Lindex <- S4Arrays:::normalize_Nindex(Nindex, seed)
    seed_dim <- dim(seed)
    nogap <- S4Arrays:::subscript_has_nogap(Lindex, seed_dim)
    if (!(is.null(dim(Rvalue)) && is.vector(Rvalue) && length(Rvalue) == 1L)) {
        selection_dim <- S4Arrays:::get_Nindex_lengths(Lindex, seed_dim)
        Rvalue <- .normarg_Rvalue(Rvalue, selection_dim)
        ## For each non-NULL subscript, keep **last** duplicate only and
        ## replace all previous duplicates with NAs.
        Lindex <- lapply(Lindex,
            function(Li) {
                if (is.null(Li))
                    return(NULL)
                Li[duplicated(Li, fromLast=TRUE)] <- NA_integer_
                Li
            })
    }
    new2("DelayedSubassign", seed=seed, Lindex=Lindex, Rvalue=Rvalue,
                             .nogap=nogap)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_noop() method
###

### Is the subassignment a no-op with respect to its "seed" slot? Note that
### even when zero array elements are being replaced, the subassignment can
### still alter the type.
setMethod("is_noop", "DelayedSubassign",
    function(x)
    {
        ## Is any array element being replaced by this subassignment?
        if (all(S4Arrays:::get_Nindex_lengths(x@Lindex, dim(x@seed)) != 0L))
            return(FALSE)
        type(x) == type(x@seed)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedSubassign

.DelayedSubassign_summary <- function(object) "Subassign"

summary.DelayedSubassign <-
    function(object, ...) .DelayedSubassign_summary(object, ...)

setMethod("summary", "DelayedSubassign", summary.DelayedSubassign)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_Mindex() and subset_DelayedSubassign()
###

### Do NOT use if 'x@Lindex' might contain duplicates! NAs are ok.
### The returned index won't contain NAs along the dimensions with no gap
### (i.e. along the dimensions for which 'x@.nogap' is TRUE).
make_Mindex <- function(index, x)
{
    stopifnot(is(x, "DelayedSubassign"),
              is.list(index),
              length(index) == length(x@Lindex))
    x_dim <- dim(x)
    lapply(seq_along(index),
        function(along) {
            i <- index[[along]]
            Li <- x@Lindex[[along]]
            if (is.null(Li))
                return(i)
            if (!is.null(i)) {
                ## match() will do the right thing if 'Li' contains NAs but
                ## NOT if it contains duplicates! This is because it will
                ## find the match to the first duplicate when we need the
                ## match to the last one.
                return(match(i, Li))
            }
            d <- x_dim[[along]]
            ## A slightly faster version of 'match(seq_len(d), Li)'. All the
            ## non-NA values in 'Li' are supposed to be >= 1 and <= d.
            m <- rep.int(NA_integer_, d)
            nonNA_idx <- which(!is.na(Li))
            m[Li[nonNA_idx]] <- seq_along(Li)[nonNA_idx]
            m
        })
}

### The returned index should never contain NAs!
.get_Lindex2_from_Mindex <- function(Mindex, nogap)
{
    lapply(seq_along(Mindex),
        function(along) {
            if (nogap[[along]])
                return(NULL)
            m <- Mindex[[along]]
            Li2 <- which(!is.na(m))
            if (length(Li2) == length(m))
                return(NULL)
            Li2
        })
}

### A more efficient version of .get_Lindex2_from_Mindex(make_Mindex(...))
### that can only be used when the right value of the subassignment is an
### ordinary vector of length 1.
### Assume that 'x@Lindex' does NOT contain NAs. Duplicates are ok.
### The returned index should never contain NAs!
.make_Lindex2 <- function(index, x)
{
    stopifnot(is(x, "DelayedSubassign"),
              is.list(index),
              length(index) == length(x@Lindex))
    lapply(seq_along(index),
        function(along) {
            if (x@.nogap[[along]])
                return(NULL)
            i <- index[[along]]
            Li <- x@Lindex[[along]]
            if (is.null(i))
                return(Li)
            Li2 <- which(i %in% Li)
            if (length(Li2) == length(i))
                return(NULL)
            Li2
        })
}

### The returned index should never contain NAs!
.get_Rindex_from_Mindex <- function(Mindex, Lindex2)
{
    lapply(seq_along(Mindex),
        function(along) {
            m <- Mindex[[along]]
            if (is.null(Lindex2[[along]]))
                return(m)
            m[!is.na(m)]
        })
}

### 'index' is assumed to be a normalized Nindex compatible with
### DelayedSubassign object 'x'.
### Return a DelayedSubassign object that represents the action of subsetting
### 'x' with 'index'. This new DelayedSubassign object is obtained by:
### - replacing 'x@Lindex' with a left index that contains strictly sorted
###   subscripts with no NAs;
### - replacing 'x@seed' with a DelayedSubset object that represents the
###   action of subsetting it with 'index';
### - if 'x@Rvalue' is an array-like object, replacing it with a DelayedSubset
###   object that represents the action of subsetting it with the index
###   returned by .get_Rindex_from_Mindex().
subset_DelayedSubassign <- function(x, index=NULL)
{
    stopifnot(is(x, "DelayedSubassign"))
    if (is.null(index))
        index <- vector("list", length=length(x@Lindex))
    ans_seed <- new2("DelayedSubset", seed=x@seed, index=index, check=FALSE)
    if (is.null(dim(x@Rvalue))) {
        ## 'x@Rvalue' is an ordinary vector (atomic or list) of length 1
        ans_Lindex <- .make_Lindex2(index, x)
        ans_Rvalue <- x@Rvalue
    } else {
        ## 'x@Rvalue' is an array-like object
        Mindex <- make_Mindex(index, x)
        ans_Lindex <- .get_Lindex2_from_Mindex(Mindex, x@.nogap)
        Rindex <- .get_Rindex_from_Mindex(Mindex, ans_Lindex)
        ans_Rvalue <- new2("DelayedSubset", seed=x@Rvalue, index=Rindex,
                                            check=FALSE)
    }
    ans_nogap <- S4Arrays:::subscript_has_nogap(ans_Lindex, dim(ans_seed))
    new2("DelayedSubassign", seed=ans_seed,
                             Lindex=ans_Lindex,
                             Rvalue=ans_Rvalue,
                             .nogap=ans_nogap,
                             check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###
### We inherit the default dim() and dimnames() methods defined for
### DelayedUnaryIsoOp derivatives, but overwite their extract_array() method.

.extract_array_from_DelayedSubassign <- function(x, index)
{
    x2 <- subset_DelayedSubassign(x, index)
    if (is.null(dim(x2@Rvalue))) {
        ## 'x2@Rvalue' is an ordinary vector (atomic or list) of length 1
        a2 <- x2@Rvalue
    } else {
        ## 'x2@Rvalue' is an array-like object
        a2 <- extract_array(x2@Rvalue@seed, x2@Rvalue@index)
    }
    if (all(x2@.nogap)) {
        if (is.null(dim(x2@Rvalue))) {
            a_dim <- S4Arrays:::get_Nindex_lengths(index, dim(x2@seed))
            a2 <- array(a2, a_dim)
        }
        return(a2)
    }
    a <- extract_array(x2@seed@seed, x2@seed@index)
    S4Arrays:::replace_by_Nindex(a, x2@Lindex, a2)
}

setMethod("extract_array", "DelayedSubassign",
    .extract_array_from_DelayedSubassign
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedSubassign",
    function(x) {
        ## We return FALSE for now.
        ## TODO: Implement this.
        FALSE
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedSubassign",
    function(x, index)
    {
        stop("NOT IMPLEMENTED YET!")
    }
)

