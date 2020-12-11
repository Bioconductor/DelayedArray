### =========================================================================
### DelayedOp objects
### -------------------------------------------------------------------------
###
### In a DelayedArray object the delayed operations are stored as a tree of
### DelayedOp objects. Each node in the tree is represented by a DelayedOp
### object. 8 types of nodes are currently supported. Each type is a concrete
### DelayedOp subclass:
###
###   Node type                        Represented operation
###   -------------------------------------------------------------------
###   DelayedOp (VIRTUAL)
###   -------------------------------------------------------------------
###   * DelayedUnaryOp (VIRTUAL)
###     o DelayedSubset                Multi-dimensional single bracket
###                                    subsetting.
###     o DelayedAperm                 Extended aperm() (can drop and/or
###                                    add ineffective dimensions).
###     o DelayedUnaryIsoOp (VIRTUAL)  Unary op that preserves the
###                                    geometry.
###       - DelayedUnaryIsoOpStack     Simple ops stacked together.
###       - DelayedUnaryIsoOpWithArgs  One op with vector-like arguments
###                                    along the dimensions of the input.
###       - DelayedSubassign           Multi-dimensional single bracket
###                                    subassignment.
###       - DelayedDimnames            Set/replace the dimnames.
###   -------------------------------------------------------------------
###   * DelayedNaryOp (VIRTUAL)
###     o DelayedNaryIsoOp             N-ary op that preserves the
###                                    geometry.
###     o DelayedAbind                 abind()
###   -------------------------------------------------------------------
###
### All the nodes are array-like objects that must comply with the "seed
### contract" i.e. they must support dim(), dimnames(), and extract_array().
###

### This virtual class and its 8 concrete subclasses are for internal use
### only and never exposed to the end user.
setClass("DelayedOp", contains="Array", representation("VIRTUAL"))

### NOT exported for now.
setGeneric("is_noop", function(x) standardGeneric("is_noop"))

### S3/S4 combo for summary.DelayedOp

.DelayedOp_summary <- function(object) sprintf("%s object", class(object))

summary.DelayedOp <- function(object, ...) .DelayedOp_summary(object, ...)

setMethod("summary", "DelayedOp", summary.DelayedOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedUnaryOp objects
###

setClass("DelayedUnaryOp",
    contains="DelayedOp",
    representation(
        "VIRTUAL",
        seed="ANY"    # The input array-like object. Expected to comply
                      # with the "seed contract".
    ),
    prototype(
        seed=new("array")
    )
)

.validate_DelayedUnaryOp <- function(x)
{
    if (length(dim(x@seed)) == 0L)
        return("the supplied seed must have dimensions")
    TRUE
}

setValidity2("DelayedUnaryOp", .validate_DelayedUnaryOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedSubset objects
###
### Delayed "Multi-dimensional single bracket subsetting".
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

### 'Nindex' must be a "multidimensional subsetting Nindex" (see
### Nindex-utils.R) or NULL.
new_DelayedSubset <- function(seed=new("array"), Nindex=NULL)
{
    index <- normalizeNindex(Nindex, seed)
    new2("DelayedSubset", seed=seed, index=index)
}

setMethod("is_noop", "DelayedSubset",
    function(x) all(S4Vectors:::sapply_isNULL(x@index))
)

### S3/S4 combo for summary.DelayedSubset

.DelayedSubset_summary <- function(object) "Subset"

summary.DelayedSubset <-
    function(object, ...) .DelayedSubset_summary(object, ...)

setMethod("summary", "DelayedSubset", summary.DelayedSubset)

### Seed contract.

setMethod("dim", "DelayedSubset",
    function(x) get_Nindex_lengths(x@index, dim(x@seed))
)

setMethod("dimnames", "DelayedSubset",
    function(x) subset_dimnames_by_Nindex(dimnames(x@seed), x@index)
)

setMethod("extract_array", "DelayedSubset",
    function(x, index)
    {
        x2 <- subset_DelayedSubset(x, index)
        extract_array(x2@seed, x2@index)
    }
)

### is_sparse() and extract_sparse_array()

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedAperm objects
###
### Delayed "Extended aperm()" (can drop and/or add ineffective dimensions).
### Note that since only "ineffective" dimensions (i.e. dimensions equal to 1)
### can be dropped or added, length is always preserved.
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

new_DelayedAperm <- function(seed=new("array"), perm=NULL)
{
    perm <- normarg_perm(perm, dim(seed))
    new2("DelayedAperm", seed=seed, perm=perm)
}

setMethod("is_noop", "DelayedAperm",
    function(x) isSequence(x@perm, length(dim(x@seed)))
)

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

### Seed contract.

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

### is_sparse() and extract_sparse_array()

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
### DelayedUnaryIsoOp objects
###
### Delayed "Unary op that preserves the geometry".
###

setClass("DelayedUnaryIsoOp",
    contains="DelayedUnaryOp",
    representation("VIRTUAL")
)

### Seed contract.
### The 3 default methods below give DelayedUnaryIsoOp derivatives a no-op
### semantic by default.
### DelayedUnaryIsoOpStack and DelayedUnaryIsoOpWithArgs objects overwrite
### this default "extract_array" method.
### DelayedDimnames objects overwrite this default "dimnames" method.
### Note that a DelayedArray object is also a DelayedUnaryIsoOp derivative
### and is considered to be the root node of the tree of DelayedOp objects
### contained in it. From a DelayedOp point of view, this root node must
### represent a no-op so DelayedArray objects inherit the 3 default methods
### below.

setMethod("dim", "DelayedUnaryIsoOp", function(x) dim(x@seed))

setMethod("dimnames", "DelayedUnaryIsoOp", function(x) dimnames(x@seed))

setMethod("extract_array", "DelayedUnaryIsoOp",
    function(x, index) extract_array(x@seed, index)
)

.set_or_check_dim <- function(x, dim)
{
    x_dim <- dim(x)
    if (is.null(x_dim)) {
        dim(x) <- dim
    } else {
        stopifnot(identical(x_dim, dim))
    }
    x
}

### is_sparse() and extract_sparse_array()
### Like the 3 default methods above (seed contract), the 2 default methods
### below also implement a no-op semantic and are also inherited by
### DelayedArray objects.

setMethod("is_sparse", "DelayedUnaryIsoOp", function(x) is_sparse(x@seed))

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedUnaryIsoOp",
    function(x, index) extract_sparse_array(x@seed, index)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedUnaryIsoOpStack objects
###
### Delayed "Unary op that preserves the geometry" where the op is made of
### simple ops stacked together.
###

setClass("DelayedUnaryIsoOpStack",
    contains="DelayedUnaryIsoOp",
    representation(
        OPS="list"  # The functions to apply to the input i.e. to the
                    # incoming array-like object. For example log
                    # or function(x) log(x + 1). It should act as an
                    # isomorphism i.e. always output an array-like
                    # object **parallel** to the input (i.e. with the
                    # same dimensions as the input).
    ),
    prototype(
        OPS=list()
    )
)

new_DelayedUnaryIsoOpStack <- function(seed=new("array"), OPS=list(),
                                       check.op=FALSE)
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    if (!is.list(OPS))
        stop(wmsg("'OPS' must be a list"))
    OPS <- lapply(OPS, match.fun)

    ans <- new2("DelayedUnaryIsoOpStack", seed=seed, OPS=OPS)
    if (check.op) {
        ## We quickly test the validity of the operation by calling type()
        ## on the returned object. This will fail if the operation cannot
        ## be applied e.g. if the user does something like:
        ##   M <- DelayedArray(matrix(character(12), ncol=3))
        ##   M2 <- log(M)
        ## The test is cheap and type() will be called anyway by show()
        ## later when the user tries to display M2. Better fail early than
        ## late!
        type(ans)  # we ignore the returned value
    }
    ans
}

### S3/S4 combo for summary.DelayedUnaryIsoOpStack

.DelayedUnaryIsoOpStack_summary <- function(object) "Unary iso op stack"

summary.DelayedUnaryIsoOpStack <-
    function(object, ...) .DelayedUnaryIsoOpStack_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOpStack", summary.DelayedUnaryIsoOpStack)

### Seed contract.
### We inherit the "dim" and "dimnames" default methods for DelayedUnaryIsoOp
### derivatives, and overwite their "extract_array" method.

setMethod("extract_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        a <- extract_array(x@seed, index)
        a_dim <- dim(a)
        for (OP in x@OPS) {
            a <- OP(a)
            ## Some operations (e.g. dnorm()) don't propagate the "dim"
            ## attribute if the input array is empty.
            a <- .set_or_check_dim(a, a_dim)
        }
        a
    }
)

### is_sparse() and extract_sparse_array()

### Make an ordinary array of the specified type and number of dimensions,
### and with a single "zero" element. The single element is the "zero"
### associated with the specified type e.g. the empty string ("") if type
### is "character", FALSE if it's "logical", etc... More generally, the
### "zero" element is whatever 'vector(type, length=1L)' produces.
.make_array_of_one_zero <- function(type, ndim)
{
    array(vector(type, length=1L), dim=rep.int(1L, ndim))
}

setMethod("is_sparse", "DelayedUnaryIsoOpStack",
    function(x)
    {
        if (!is_sparse(x@seed))
            return(FALSE)
        ## Structural sparsity will be propagated if the operations in
        ## x@OPS preserve the zeroes. To find out whether zeroes are preserved
        ## or not, we replace the current seed with an array of one "zero",
        ## that is, with an ordinary array of the same number of dimensions
        ## and type as the seed, but with a single "zero" element. Then we
        ## apply the operations in x@OPS to it and see whether the zero was
        ## preserved or not.
        seed_ndim <- length(dim(x@seed))
        x@seed <- .make_array_of_one_zero(type(x@seed), seed_ndim)
        a0 <- extract_array(x, rep.int(list(1L), seed_ndim))
        as.vector(a0) == vector(type(a0), length=1L)
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        ## Assuming that the caller respected "extract_sparse_array() Terms
        ## of Use" (see SparseArraySeed-class.R), 'is_sparse(x)' should be
        ## TRUE so we can assume that the operations in x@OPS preserve the
        ## zeroes and thus only need to apply them to the nonzero data.
        sas <- extract_sparse_array(x@seed, index)
        sas_nzdata <- sas@nzdata
        for (OP in x@OPS)
            sas_nzdata <- OP(sas_nzdata)
        sas@nzdata <- sas_nzdata
        sas
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedUnaryIsoOpWithArgs objects
###
### Delayed "Unary op with arguments that preserves the geometry".
### Here the op can have vector-like arguments along the dimensions of the
### input.
###

setClass("DelayedUnaryIsoOpWithArgs",
    contains="DelayedUnaryIsoOp",
    representation(
        OP="function",     # The function to apply to the input i.e. to the
                           # incoming array-like object. For example `+` or
                           # log. It should act as an isomorphism i.e. always
                           # output an array-like object **parallel** to the
                           # input (i.e. with the same dimensions as the input).

        Largs="list",      # Left arguments to OP i.e. arguments to place
                           # before the input array in the function call.
        Rargs="list",      # Right arguments to OP i.e. arguments to place
                           # after the input array in the function call.

        Lalong="integer",  # One integer (or NA) per left argument indicating
                           # which dimension of the input array the argument
                           # is parallel to.
        Ralong="integer"   # One integer (or NA) per right argument indicating
                           # which dimension of the input array the argument
                           # is parallel to.
    ),
    prototype(
        OP=identity
    )
)

.normarg_Lalong_or_Ralong <- function(Lalong, Largs, seed_dim)
{
    if (identical(Lalong, NA))
        return(rep.int(NA_integer_, length(Largs)))
    if (!(is.numeric(Lalong) && length(Lalong) == length(Largs)))
        stop(wmsg("'Lalong' and 'Ralong' must be integer vectors ",
                  "parallel to 'Largs' and 'Rargs', respectively"))
    if (!is.integer(Lalong))
        Lalong <- as.integer(Lalong)
    nonNA_idx <- which(!is.na(Lalong))
    nonNA_Lalong <- Lalong[nonNA_idx]
    if (S4Vectors:::anyMissingOrOutside(nonNA_Lalong, 1L, length(seed_dim)))
        stop(wmsg("all non-NA values in 'Lalong' and 'Ralong' must ",
                  "be >= 1 and <= 'length(dim(seed))'"))
    if (any(Lalong != 1L, na.rm=TRUE))
        stop(wmsg("arguments in 'Largs' and 'Rargs' can only go along ",
                  "with the first dimension at the moment"))
    ok <- elementNROWS(Largs[nonNA_idx]) == seed_dim[nonNA_Lalong]
    if (!all(ok))
        stop(wmsg("some arguments in 'Largs' and/or 'Rargs' are not ",
                  "parallel to the dimension that they go along with"))
    Lalong
}

new_DelayedUnaryIsoOpWithArgs <- function(seed=new("array"),
                                          OP=identity,
                                          Largs=list(), Rargs=list(),
                                          Lalong=NA, Ralong=NA,
                                          check.op=FALSE)
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    stopifnot(is.list(Largs), is.list(Rargs))
    Lalong <- .normarg_Lalong_or_Ralong(Lalong, Largs, seed_dim)
    Ralong <- .normarg_Lalong_or_Ralong(Ralong, Rargs, seed_dim)

    OP <- match.fun(OP)

    ans <- new2("DelayedUnaryIsoOpWithArgs", seed=seed,
                                             OP=OP,
                                             Largs=Largs, Rargs=Rargs,
                                             Lalong=Lalong, Ralong=Ralong)
    if (check.op)
        type(ans)  # we ignore the returned value
    ans
}

### S3/S4 combo for summary.DelayedUnaryIsoOpWithArgs

.DelayedUnaryIsoOpWithArgs_summary <- function(object) "Unary iso op with args"

summary.DelayedUnaryIsoOpWithArgs <-
    function(object, ...) .DelayedUnaryIsoOpWithArgs_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOpWithArgs",
    summary.DelayedUnaryIsoOpWithArgs
)

### Seed contract.
### We inherit the "dim" and "dimnames" default methods for DelayedUnaryIsoOp
### derivatives, and overwite their "extract_array" method.

subset_args <- function(args, along, index)
{
    subset_arg <- function(arg, MARGIN) {
        if (is.na(MARGIN))
            return(arg)
        i <- index[[MARGIN]]
        if (is.null(i))
            return(arg)
        extractROWS(arg, i)
    }
    mapply(subset_arg, args, along, SIMPLIFY=FALSE, USE.NAMES=FALSE)
}

setMethod("extract_array", "DelayedUnaryIsoOpWithArgs",
    function(x, index)
    {
        a <- extract_array(x@seed, index)

        ## Subset the left and right arguments that go along with a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        ans <- do.call(x@OP, c(Largs, list(a), Rargs))

        ## Some operations (e.g. dnorm()) don't propagate the "dim" attribute
        ## if the input array is empty.
        .set_or_check_dim(ans, dim(a))
    }
)

### is_sparse() and extract_sparse_array()

### DelayedUnaryIsoOpWithArgs objects are NOT considered to propagate
### structural sparsity at the moment. However it would be nice if
### things like 'M * runif(nrow(M))' or 'M / runif(nrow(M))' propagated
### sparsity. These are simplified versions of the following use case by
### Aaron:
###   library(TENxBrainData)
###   sce <- TENxBrainData()
###   sf <- runif(ncol(sce))
###   lcounts <- log2(t(t(counts(sce))/sf) + 1)
### 'lcounts' should be considered sparse but right now it's not!
### TODO: The "is_sparse" method below should be able to propagate sparsity
### of 'A * v', 'v * A', and 'A / v',  when 'length(v)' is 'nrow(A)' and
### 'v' does not contain infinite or NA or NaN values (in the multiplication
### case) and no zero or NA or NaN values (in the division case).
setMethod("is_sparse", "DelayedUnaryIsoOpWithArgs", function(x) FALSE)

setMethod("extract_sparse_array", "DelayedUnaryIsoOpWithArgs",
    function(x, index)
        stop(wmsg("extract_sparse_array() is not supported ",
                  "on DelayedUnaryIsoOpWithArgs objects"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedSubassign objects
###
### Delayed "Multi-dimensional single bracket subassignment".
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
### Nindex-utils.R) or NULL.
new_DelayedSubassign <- function(seed=new("array"), Nindex=NULL, Rvalue=NA)
{
    Lindex <- normalizeNindex(Nindex, seed)
    seed_dim <- dim(seed)
    nogap <- subscript_has_nogap(Lindex, seed_dim)
    if (!(is.null(dim(Rvalue)) && is.vector(Rvalue) && length(Rvalue) == 1L)) {
        selection_dim <- get_Nindex_lengths(Lindex, seed_dim)
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

### Is the subassignment a no-op with respect to its "seed" slot? Note that
### even when zero array elements are being replaced, the subassignment can
### still alter the type.
setMethod("is_noop", "DelayedSubassign",
    function(x)
    {
        ## Is any array element being replaced by this subassignment?
        if (all(get_Nindex_lengths(x@Lindex, dim(x@seed)) != 0L))
            return(FALSE)
        type(x) == type(x@seed)
    }
)

### S3/S4 combo for summary.DelayedSubassign

.DelayedSubassign_summary <- function(object) "Subassign"

summary.DelayedSubassign <-
    function(object, ...) .DelayedSubassign_summary(object, ...)

setMethod("summary", "DelayedSubassign", summary.DelayedSubassign)

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
    ans_nogap <- subscript_has_nogap(ans_Lindex, dim(ans_seed))
    new2("DelayedSubassign", seed=ans_seed,
                             Lindex=ans_Lindex,
                             Rvalue=ans_Rvalue,
                             .nogap=ans_nogap,
                             check=FALSE)
}

### Seed contract.
### We inherit the "dim" and "dimnames" default methods for DelayedUnaryIsoOp
### derivatives, and overwite their "extract_array" method.

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
            a_dim <- get_Nindex_lengths(index, dim(x2@seed))
            a2 <- array(a2, a_dim)
        }
        return(a2)
    }
    a <- extract_array(x2@seed@seed, x2@seed@index)
    replace_by_Nindex(a, x2@Lindex, a2)
}

setMethod("extract_array", "DelayedSubassign",
    .extract_array_from_DelayedSubassign
)

### is_sparse() and extract_sparse_array()

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedDimnames objects
###
### Delayed "Set dimnames".
###

### Used in unit tests!
.INHERIT_FROM_SEED <- -1L

setClass("DelayedDimnames",
    contains="DelayedUnaryIsoOp",
    representation(
        dimnames="list"  # List with one list element per dimension in
                         # the input. Each list element must be NULL,
                         # or a character vector, or special value
                         # .INHERIT_FROM_SEED
    ),
    prototype(
        dimnames=list(.INHERIT_FROM_SEED)
    )
)

.validate_DelayedDimnames <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)

    ## 'dimnames' slot.
    if (length(x@dimnames) != seed_ndim)
        return(paste0("'x@dimnames' must have one list element per ",
                      "dimension in 'x@seed'"))
    ok <- mapply(function(dn, d) {
                     identical(dn, .INHERIT_FROM_SEED) ||
                     is.null(dn) ||
                     is.character(dn) && length(dn) == d
                 },
                 x@dimnames, seed_dim,
                 SIMPLIFY=FALSE, USE.NAMES=FALSE)
    if (!all(unlist(ok)))
        return(paste0("each list element in 'x@dimnames' must be NULL, ",
                      "or a character vector of length the extent of ",
                      "the corresponding dimension, or special value ",
                      .INHERIT_FROM_SEED))
    TRUE
}

setValidity2("DelayedDimnames", .validate_DelayedDimnames)

new_DelayedDimnames <- function(seed=new("array"), dimnames=.INHERIT_FROM_SEED)
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    if (identical(dimnames, .INHERIT_FROM_SEED)) {
        dimnames <- rep.int(list(.INHERIT_FROM_SEED), seed_ndim)
    } else {
        dimnames <- normarg_dimnames(dimnames, seed_dim)
        seed_dimnames <- dimnames(seed)
        dimnames <- lapply(setNames(seq_len(seed_ndim), names(dimnames)),
                           function(along) {
                               dn <- dimnames[[along]]
                               seed_dn <- seed_dimnames[[along]]
                               ## Let's play nice with seeds that return
                               ## dimnames() that are not NULL or a character
                               ## vector e.g.
                               ##   library(GDSArray)
                               ##   gds <- seqExampleFileName("gds")
                               ##   ga1 <- GDSArray(gds, "genotype/data")
                               ##   sapply(dimnames(ga1), class)
                               ##   #  variant.id   sample.id   ploidy.id
                               ##   # "character" "character"   "integer"
                               if (!(is.null(seed_dn) || is.character(seed_dn)))
                                   seed_dn <- as.character(seed_dn)
                               if (identical(dn, seed_dn))
                                   return(.INHERIT_FROM_SEED)
                               dn
                           })
    }
    new2("DelayedDimnames", seed=seed, dimnames=dimnames)
}

setMethod("is_noop", "DelayedDimnames",
    function(x)
    {
        ok <- vapply(x@dimnames, identical, logical(1), .INHERIT_FROM_SEED,
                     USE.NAMES=FALSE)
        all(ok) && identical(names(x@dimnames), names(dimnames(x@seed)))
    }
)

### S3/S4 combo for summary.DelayedDimnames

.DelayedDimnames_summary <- function(object) "Set dimnames"

summary.DelayedDimnames <-
    function(object, ...) .DelayedDimnames_summary(object, ...)

setMethod("summary", "DelayedDimnames", summary.DelayedDimnames)

### Seed contract.
### We inherit the "dim" and "extract_array" default methods for
### DelayedUnaryIsoOp derivatives, and overwite their "dimnames" method.

.get_DelayedDimnames_dimnames <- function(x)
{
    x_dimnames <- x@dimnames
    seed_dimnames <- dimnames(x@seed)
    ans <- lapply(setNames(seq_along(x_dimnames), names(x_dimnames)),
                  function(along) {
                      dn <- x_dimnames[[along]]
                      if (identical(dn, .INHERIT_FROM_SEED))
                          dn <- seed_dimnames[[along]]
                      dn
                  })
    simplify_NULL_dimnames(ans)
}

setMethod("dimnames", "DelayedDimnames", .get_DelayedDimnames_dimnames)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedNaryOp objects
###

setClass("DelayedNaryOp",
    contains="DelayedOp",
    representation(
        "VIRTUAL",
        seeds="list"  # The input array-like objects. Each object is
                      # expected to comply with the "seed contract".
    ),
    prototype(
        seeds=list(new("array"))
    )
)

.validate_DelayedNaryOp <- function(x)
{
    if (length(x@seeds) == 0L)
        return("'x@seeds' cannot be empty")
    TRUE
}

setValidity2("DelayedNaryOp", .validate_DelayedNaryOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedNaryIsoOp objects
###
### Delayed "N-ary op that preserves the geometry".
### The input objects must be "conformable" array-like objects i.e. they all
### must have the same dimensions.
###

setClass("DelayedNaryIsoOp",
    contains="DelayedNaryOp",
    representation(
        OP="function",  # The function to use to combine the input objects.
                        # Should act as an isomorphism i.e. always return an
                        # array-like object **parallel** to the input objects
                        # (i.e. with the same dimensions).

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        OP=identity
    )
)

.arrays_are_conformable <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    first_ndim <- ndims[[1L]]
    if (!all(ndims == first_ndim))
        return(FALSE)
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return(FALSE)
    dims <- matrix(tmp, ncol=length(objects))
    first_dim <- dims[ , 1L]
    all(dims == first_dim)
}

.validate_DelayedNaryIsoOp <- function(x)
{
    ## 'seeds' slot.
    if (!.arrays_are_conformable(x@seeds))
        return("'x@seeds' must be a list of conformable array-like objects")
    TRUE
}

setValidity2("DelayedNaryIsoOp", .validate_DelayedNaryIsoOp)

new_DelayedNaryIsoOp <- function(OP=identity, seed=new("array"), ...,
                                 Rargs=list())
{
    OP <- match.fun(OP)
    seeds <- unname(list(seed, ...))
    if (!.arrays_are_conformable(seeds))
        stop(wmsg("non-conformable array-like objects"))
    new2("DelayedNaryIsoOp", seeds=seeds, OP=OP, Rargs=Rargs, check=FALSE)
}

### S3/S4 combo for summary.DelayedNaryIsoOp

.DelayedNaryIsoOp_summary <- function(object) "N-ary iso op"

summary.DelayedNaryIsoOp <-
    function(object, ...) .DelayedNaryIsoOp_summary(object, ...)

setMethod("summary", "DelayedNaryIsoOp", summary.DelayedNaryIsoOp)

### Seed contract.

setMethod("dim", "DelayedNaryIsoOp", function(x) dim(x@seeds[[1L]]))

setMethod("dimnames", "DelayedNaryIsoOp",
    function(x) get_first_non_NULL_dimnames(x@seeds)
)

setMethod("extract_array", "DelayedNaryIsoOp",
    function(x, index)
    {
        arrays <- lapply(x@seeds, extract_array, index)
        do.call(x@OP, c(arrays, x@Rargs))
    }
)

### is_sparse() and extract_sparse_array()

setMethod("is_sparse", "DelayedNaryIsoOp",
    function(x)
    {
        ok <- vapply(x@seeds, is_sparse, logical(1), USE.NAMES=FALSE)
        if (!all(ok))
            return(FALSE)
        if (length(x@Rargs) != 0L)
            return(FALSE)
        ## Structural sparsity will be propagated if the operation in
        ## x@OP preserves the zeroes. To find out whether zeroes are preserved
        ## or not, we replace each current seed with an array of one "zero",
        ## that is, with an ordinary array of the same number of dimensions
        ## and type as the seed, but with a single "zero" element. Then we
        ## apply the n-ary operation in x@OP to them and see whether the
        ## zero were preserved or not.
        seed_ndim <- length(dim(x@seeds[[1L]]))
        x@seeds <- lapply(x@seeds,
            function(seed) .make_array_of_one_zero(type(seed), seed_ndim))
        a0 <- extract_array(x, rep.int(list(1L), seed_ndim))
        as.vector(a0) == vector(type(a0), length=1L)
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedNaryIsoOp",
    function(x, index)
    {
        stop("NOT IMPLEMENTED YET!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedAbind objects
###
### Delayed "abind()".
###

setClass("DelayedAbind",
    contains="DelayedNaryOp",
    representation(
        along="integer"  # Single integer indicating the dimension along
                         # which to bind the input objects.
    ),
    prototype(
        along=1L
    )
)

.validate_DelayedAbind <- function(x)
{
    ## 'along' slot.
    if (!(isSingleInteger(x@along) && x@along >= 1L))
        return("'x@along' must be a single positive integer")
    ndim <- length(dim(x@seeds[[1L]]))
    if (ndim < x@along)
        return(paste0("the array-like objects to bind must have at least ",
                      x@along, " dimensions for this binding operation"))

    dims <- get_dims_to_bind(x@seeds, x@along)
    if (is.character(dims))
        return(dims)
    TRUE
}

setValidity2("DelayedAbind", .validate_DelayedAbind)

new_DelayedAbind <- function(seeds, along)
{
    new2("DelayedAbind", seeds=seeds, along=along)
}

setMethod("is_noop", "DelayedAbind",
    function(x) length(x@seeds) == 1L
)

### S3/S4 combo for summary.DelayedAbind

.DelayedAbind_summary <-
    function(object) sprintf("Abind (along=%d)", object@along)

summary.DelayedAbind <-
    function(object, ...) .DelayedAbind_summary(object, ...)

setMethod("summary", "DelayedAbind", summary.DelayedAbind)

### Seed contract.

.get_DelayedAbind_dim <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dims_along(dims, x@along)
}

setMethod("dim", "DelayedAbind", .get_DelayedAbind_dim)

.get_DelayedAbind_dimnames <- function(x)
{
    dims <- get_dims_to_bind(x@seeds, x@along)
    combine_dimnames_along(x@seeds, dims, x@along)
}

setMethod("dimnames", "DelayedAbind", .get_DelayedAbind_dimnames)

.extract_array_from_DelayedAbind <- function(x, index)
{
    i <- index[[x@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        tmp <- lapply(x@seeds, extract_array, index)
        ## Bind the ordinary arrays in 'tmp'.
        ans <- do.call(simple_abind, c(tmp, list(along=x@along)))
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- get_part_index(i, breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[x@along]] <- split_part_idx[[s]]
        extract_array(x@seeds[[s]], index)
    }
    tmp <- lapply(seq_along(x@seeds), FUN)

    ## Bind the ordinary arrays in 'tmp'.
    ans <- do.call(simple_abind, c(tmp, list(along=x@along)))

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector("list", length=length(index))
    Nindex[[x@along]] <- get_rev_index(part_idx)
    extract_array(ans, Nindex)
}

setMethod("extract_array", "DelayedAbind", .extract_array_from_DelayedAbind)

### is_sparse() and extract_sparse_array()

setMethod("is_sparse", "DelayedAbind",
    function(x)
    {
        all(vapply(x@seeds, is_sparse, logical(1), USE.NAMES=FALSE))
    }
)

.extract_sparse_array_DelayedAbind <- function(x, index)
{
    i <- index[[x@along]]

    if (is.null(i)) {
        ## This is the easy situation.
        tmp <- lapply(x@seeds, extract_sparse_array, index)
        ## Bind the ordinary arrays in 'tmp'.
        ans <- abind_SparseArraySeed_objects(tmp, x@along)
        return(ans)
    }

    ## From now on 'i' is a vector of positive integers.
    dims <- get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- get_part_index(i, breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    FUN <- function(s) {
        index[[x@along]] <- split_part_idx[[s]]
        extract_sparse_array(x@seeds[[s]], index)
    }
    tmp <- lapply(seq_along(x@seeds), FUN)

    ## Bind the ordinary arrays in 'tmp'.
    ans <- abind_SparseArraySeed_objects(tmp, x@along)

    ## Reorder the rows or columns in 'ans'.
    Nindex <- vector("list", length=length(index))
    Nindex[[x@along]] <- get_rev_index(part_idx)
    extract_sparse_array(ans, Nindex)
}

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedAbind",
    .extract_sparse_array_DelayedAbind
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### In DelayedArray 0.5.24, the SeedDimPicker, ConformableSeedCombiner, and
### SeedBinder classes were renamed DelayedAperm, DelayedNaryIsoOp, and
### DelayedAbind, respectively.
### DelayedArray objects serialized with DelayedArray < 0.5.24 might contain
### instances of these old classes nested in their "seed" slot so we need to
### keep them around for now.
###

setClass("SeedDimPicker", contains="DelayedAperm")
setClass("ConformableSeedCombiner", contains="DelayedNaryIsoOp")
setClass("SeedBinder", contains="DelayedAbind")

setMethod("updateObject", "DelayedOp",
    function(object, ..., verbose=FALSE)
    {
        if (.hasSlot(object, "seed")) {
            object@seed <- updateObject(object@seed, ..., verbose=verbose)
        }
        if (.hasSlot(object, "seeds")) {
            object@seeds <- lapply(object@seeds,
                function(seed) updateObject(seed, ..., verbose=verbose))
        }
        object
    }
)

setMethod("updateObject", "SeedDimPicker",
    function(object, ..., verbose=FALSE)
    {
        object <- new2("DelayedAperm", seed=object@seed,
                                       perm=object@dim_combination)
        callNextMethod()
    }
)

setMethod("updateObject", "ConformableSeedCombiner",
    function(object, ..., verbose=FALSE)
    {
        object <- new2("DelayedNaryIsoOp", seeds=object@seeds,
                                               OP=object@COMBINING_OP,
                                               Rargs=object@Rargs)
        callNextMethod()
    }
)

setMethod("updateObject", "SeedBinder",
    function(object, ..., verbose=FALSE)
    {
        class(object) <- "DelayedAbind"
        callNextMethod()
    }
)

