### =========================================================================
### DelayedOp objects
### -------------------------------------------------------------------------
###
### In a DelayedArray object the delayed operations are stored as a tree of
### DelayedOp objects. Each node in the tree is represented by a DelayedOp
### object. 6 types of nodes are currently supported. Each type is a concrete
### DelayedOp subclass:
###
###   Node type    Outdegree  Operation
###   -------------------------------------------------------------------
###   DelayedSubset        1  Multi-dimensional single bracket subsetting
###   DelayedAperm         1  Extended aperm() (can drop dimensions)
###   DelayedUnaryIsoOp    1  Unary op that preserves the geometry
###   DelayedDimnames      1  Set dimnames
###   DelayedNaryIsoOp     N  N-ary op that preserves the geometry
###   DelayedAbind         N  abind()
###
### All the nodes are array-like objects that must comply with the "seed
### contract" i.e. they must support dim(), dimnames(), and extract_array().
###

### This virtual class and its 6 concrete subclasses are for internal use
### only and are not exported.
setClass("DelayedOp", contains="Array", representation("VIRTUAL"))

### NOT exported for now.
setGeneric("isNoOp", function(x) standardGeneric("isNoOp"))

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
        return(wmsg2("the supplied seed must have dimensions"))
    TRUE
}

setValidity2("DelayedUnaryOp", .validate_DelayedUnaryOp)

### Seed contract.
### Each DelayedUnaryOp derivative inherits these 3 default methods and will
### typically overwrite at least one of them (otherwise it would be a no-op).

setMethod("dim", "DelayedUnaryOp", function(x) dim(x@seed))

setMethod("dimnames", "DelayedUnaryOp", function(x) dimnames(x@seed))

setMethod("extract_array", "DelayedUnaryOp",
    function(x, index) extract_array(x@seed, index)
)


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
        return(wmsg2("'x@seeds' cannot be empty"))
    TRUE
}

setValidity2("DelayedNaryOp", .validate_DelayedNaryOp)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedSubset objects
###
### Delayed "Multi-dimensional single bracket subsetting".
###

setClass("DelayedSubset",
    contains="DelayedUnaryOp",
    representation(
        index="list"  # List of subscripts as positive integer vectors, one
                      # per dimension in the input. *Missing* list elements
                      # are allowed and represented by NULLs.
    ),
    prototype(
        index=list(NULL)
    )
)

.validate_DelayedSubset <- function(x)
{
    ## 'index' slot.
    if (length(x@index) != length(dim(x@seed)))
        return(wmsg2("'x@index' must have one list element per dimension ",
                     "in 'x@seed'"))
    if (!is.null(names(x@index)))
        return(wmsg2("'x@index' should not have names"))
    ok <- lapply(x@index,
              function(i) {is.null(i) || is.integer(i) && is.null(names(i))})
    if (!all(unlist(ok)))
        return(wmsg2("each list element in 'x@index' must be NULL ",
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

### 'Nindex' must be a "multidimensional subsetting Nindex" (see utils.R)
### or NULL.
new_DelayedSubset <- function(seed=new("array"), Nindex=NULL)
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    if (is.null(Nindex)) {
        index <- rep.int(list(NULL), seed_ndim)
    } else {
        stopifnot(is.list(Nindex), length(Nindex) == seed_ndim)
        ## Normalize 'Nindex' i.e. check and turn its non-NULL list elements
        ## into positive integer vectors.
        seed_dimnames <- dimnames(seed)
        index <- lapply(seq_len(seed_ndim),
            function(along) {
                subscript <- Nindex[[along]]
                if (is.null(subscript))
                    return(NULL)
                d <- seed_dim[[along]]
                i <- normalizeSingleBracketSubscript2(subscript,
                                    d, seed_dimnames[[along]])
                if (isSequence(i, of.length=d))
                    return(NULL)
                i
            })
    }
    new2("DelayedSubset", seed=seed, index=index)
}

setMethod("isNoOp", "DelayedSubset",
    function(x) all(S4Vectors:::sapply_isNULL(x@index))
)

### S3/S4 combo for summary.DelayedSubset

.DelayedSubset_summary <- function(object) "Subset"

summary.DelayedSubset <-
    function(object, ...) .DelayedSubset_summary(object, ...)

setMethod("summary", "DelayedSubset", summary.DelayedSubset)

### Seed contract.

.get_DelayedSubset_dim <- function(x) get_Nindex_lengths(x@index, dim(x@seed))

setMethod("dim", "DelayedSubset", .get_DelayedSubset_dim)

setMethod("dimnames", "DelayedSubset",
    function(x) subset_dimnames_by_Nindex(dimnames(x@seed), x@index)
)

.extract_array_from_DelayedSubset <- function(x, index)
{
    x <- subset_DelayedSubset(x, index)
    extract_array(x@seed, x@index)
}

setMethod("extract_array", "DelayedSubset", .extract_array_from_DelayedSubset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedAperm objects
###
### Delayed "Extended aperm()" (can drop dimensions).
### Note that only "ineffective" dimensions can be dropped (i.e. dimensions
### equal to 1, so dropping them preserves the length).
###

setClass("DelayedAperm",
    contains="DelayedUnaryOp",
    representation(
        perm="integer"  # Index into dim(seed) describing the *rearrangement*
                        # of the dimensions i.e. which dimensions of the input
                        # to keep and in which order.
    ),
    prototype(
        perm=1L
    )
)

### The seed is referred to as 'a' in the error messages. This makes them
### more meaningful for the end user in the context of calling aperm().
.validate_DelayedAperm <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)

    ## 'perm' slot.
    if (length(x@perm) == 0L)
        return(wmsg2("'perm' cannot be an empty vector"))
    if (S4Vectors:::anyMissingOrOutside(x@perm, 1L, seed_ndim))
        return(wmsg2("all values in 'perm' must be >= 1 ",
                     "and <= 'length(dim(a))'"))
    if (anyDuplicated(x@perm))
        return(wmsg2("'perm' cannot have duplicates"))
    if (!all(seed_dim[-x@perm] == 1L))
        return(wmsg2("only dimensions equal to 1 can be dropped"))
    TRUE
}

setValidity2("DelayedAperm", .validate_DelayedAperm)

new_DelayedAperm <- function(seed, perm=NULL)
{
    perm <- normarg_perm(perm, dim(seed))
    new2("DelayedAperm", seed=seed, perm=perm)
}

setMethod("isNoOp", "DelayedAperm",
    function(x) identical(x@perm, seq_along(dim(x@seed)))
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
    seed_dim[x@perm]
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

.extract_array_from_DelayedAperm <- function(x, index)
{
    seed_dim <- dim(x@seed)
    seed_index <- rep.int(list(1L), length(seed_dim))
    seed_index[x@perm] <- index
    a <- extract_array(x@seed, seed_index)
    dim(a) <- dim(a)[sort(x@perm)]
    aperm(a, perm=rank(x@perm))
}

setMethod("extract_array", "DelayedAperm",
    .extract_array_from_DelayedAperm
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedUnaryIsoOp objects
###
### Delayed "Unary op that preserves the geometry".
###

setClass("DelayedUnaryIsoOp",
    contains="DelayedUnaryOp",
    representation(
        OP="function",     # The function to apply to the input (e.g. `+` or
                           # log). It should act as an isomorphism i.e. always
                           # return an array-like object *parallel* to the
                           # input (i.e. with the same dimensions as the input).

        Largs="list",      # Left arguments to OP i.e. arguments to place
                           # before the input array in the function call.
        Rargs="list",      # Right arguments to OP i.e. arguments to place
                           # after the input array in the function call.

        Lalong="integer",  # One integer (or NA) per left argument indicating
                           # what dimension the argument is parallel to.
        Ralong="integer"   # One integer (or NA) per right argument indicating
                           # what dimension the argument is parallel to.
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
        stop(wmsg("'Lalong' must be an integer vector parallel to 'Largs'"))
    if (!is.integer(Lalong))
        Lalong <- as.integer(Lalong)
    non_na_idx <- which(!is.na(Lalong))
    non_na_Lalong <- Lalong[non_na_idx]
    if (S4Vectors:::anyMissingOrOutside(non_na_Lalong, 1L, length(seed_dim)))
        stop(wmsg("all non-NA values in 'Lalong' and 'Ralong' must ",
                  "be >= 1 and <= 'length(dim(seed))'"))
    if (any(Lalong != 1L, na.rm=TRUE))
        stop(wmsg("arguments in 'Largs' and 'Rargs' can only go along ",
                  "first dimension at the moment"))
    ok <- elementNROWS(Largs[non_na_idx]) == seed_dim[non_na_Lalong]
    if (!all(ok))
        stop(wmsg("some arguments in 'Largs' and/or 'Rargs' are not ",
                  "parallel to the dimension that they go along"))
    Lalong
}

new_DelayedUnaryIsoOp <- function(seed=new("array"),
                                  OP=identity, Largs=list(), Rargs=list(),
                                  Lalong=NA, Ralong=NA)
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    stopifnot(is.list(Largs), is.list(Rargs))
    Lalong <- .normarg_Lalong_or_Ralong(Lalong, Largs, seed_dim)
    Ralong <- .normarg_Lalong_or_Ralong(Ralong, Rargs, seed_dim)

    OP <- match.fun(OP)

    new2("DelayedUnaryIsoOp", seed=seed, OP=OP, Largs=Largs, Rargs=Rargs,
                              Lalong=Lalong, Ralong=Ralong)
}

### S3/S4 combo for summary.DelayedUnaryIsoOp

.DelayedUnaryIsoOp_summary <- function(object) "Unary iso op"

summary.DelayedUnaryIsoOp <-
    function(object, ...) .DelayedUnaryIsoOp_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOp", summary.DelayedUnaryIsoOp)

### Seed contract.

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

setMethod("extract_array", "DelayedUnaryIsoOp",
    function(x, index)
    {
        a <- extract_array(x@seed, index)

        ## Subset the left and right arguments that go along a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        ans <- do.call(x@OP, c(Largs, list(a), Rargs))

        ## Some operations (e.g. dnorm()) don't propagate the "dim" attribute
        ## if the input array is empty.
        a_dim <- dim(a)
        ans_dim <- dim(ans)
        if (is.null(ans_dim)) {
            dim(ans) <- a_dim  # propagate the "dim" attribute
        } else {
            stopifnot(identical(ans_dim, a_dim))  # sanity check
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedDimnames objects
###
### Delayed "Set dimnames".
###

.INHERIT_FROM_SEED <- -1L

setClass("DelayedDimnames",
    contains="DelayedUnaryOp",
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
        return(wmsg2("'x@dimnames' must have one list element per dimension ",
                     "in 'x@seed'"))
    ok <- mapply(function(dn, d) {
                     identical(dn, .INHERIT_FROM_SEED) ||
                     is.null(dn) ||
                     is.character(dn) && length(dn) == d
                 },
                 x@dimnames, seed_dim,
                 SIMPLIFY=FALSE, USE.NAMES=FALSE)
    if (!all(unlist(ok)))
        return(wmsg2("each list element in 'x@dimnames' must be NULL, ",
                     "or a character vector of length the extent of ",
                     "the corresponding dimension, or special value ",
                     .INHERIT_FROM_SEED))
    TRUE
}

setValidity2("DelayedDimnames", .validate_DelayedDimnames)

### TODO: Also make sure that each 'dimnames' list element is either NULL or
### a character vector of the correct length.
.normalize_dimnames <- function(dimnames, ndim)
{
    if (is.null(dimnames))
        return(vector("list", length=ndim))
    if (!is.list(dimnames))
        stop("the supplied dimnames must be a list")
    if (length(dimnames) > ndim)
        stop(wmsg("the supplied dimnames is longer ",
                  "than the number of dimensions"))
    if (length(dimnames) < ndim)
        length(dimnames) <- ndim
    dimnames
}

new_DelayedDimnames <- function(seed=new("array"),
                                dimnames=list(.INHERIT_FROM_SEED))
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    dimnames <- .normalize_dimnames(dimnames, seed_ndim)
    seed_dimnames <- dimnames(seed)
    dimnames <- lapply(seq_len(seed_ndim),
                       function(along) {
                           dn <- dimnames[[along]]
                           if (identical(dn, seed_dimnames[[along]]))
                               return(.INHERIT_FROM_SEED)
                           dn
                       })
    new2("DelayedDimnames", seed=seed, dimnames=dimnames)
}

setMethod("isNoOp", "DelayedDimnames",
    function(x)
        all(vapply(x@dimnames, identical, logical(1), .INHERIT_FROM_SEED))
)

### S3/S4 combo for summary.DelayedDimnames

.DelayedDimnames_summary <- function(object) "Set dimnames"

summary.DelayedDimnames <-
    function(object, ...) .DelayedDimnames_summary(object, ...)

setMethod("summary", "DelayedDimnames", summary.DelayedDimnames)

### Seed contract.

.get_DelayedDimnames_dimnames <- function(x)
{
    x_dimnames <- x@dimnames
    seed_dimnames <- dimnames(x@seed)
    ans <- lapply(seq_along(x_dimnames),
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
                        # array-like object *parallel* to the input objects
                        # (i.e. with the same dimensions).

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        OP=identity
    )
)

.objects_are_conformable_arrays <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    first_ndim <- ndims[[1L]]
    if (!all(ndims == first_ndim))
        return(FALSE)
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return(FALSE)
    dims <- matrix(tmp, nrow=first_ndim)
    first_dim <- dims[ , 1L]
    all(dims == first_dim)
}

.validate_DelayedNaryIsoOp <- function(x)
{
    ## 'seeds' slot.
    if (!.objects_are_conformable_arrays(x@seeds))
        return(wmsg2("'x@seeds' must be a list of conformable ",
                     "array-like objects"))
    TRUE
}

setValidity2("DelayedNaryIsoOp", .validate_DelayedNaryIsoOp)

new_DelayedNaryIsoOp <- function(seed=new("array"), ...,
                                 OP=identity, Rargs=list())
{
    seeds <- unname(list(seed, ...))
    OP <- match.fun(OP)
    new2("DelayedNaryIsoOp", seeds=seeds, OP=OP, Rargs=Rargs)
}

### S3/S4 combo for summary.DelayedNaryIsoOp

.DelayedNaryIsoOp_summary <- function(object) "N-ary iso op"

summary.DelayedNaryIsoOp <-
    function(object, ...) .DelayedNaryIsoOp_summary(object, ...)

setMethod("summary", "DelayedNaryIsoOp", summary.DelayedNaryIsoOp)

### Seed contract.

setMethod("dim", "DelayedNaryIsoOp", function(x) dim(x@seeds[[1L]]))

setMethod("dimnames", "DelayedNaryIsoOp",
    function(x) simplify_NULL_dimnames(combine_dimnames(x@seeds))
)

setMethod("extract_array", "DelayedNaryIsoOp",
    function(x, index)
    {
        arrays <- lapply(x@seeds, extract_array, index)
        do.call(x@OP, c(arrays, x@Rargs))
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
        return(wmsg2("'x@along' must be a single positive integer"))
    ndim <- length(dim(x@seeds[[1L]]))
    if (ndim < x@along)
        return(wmsg2("the array-like objects to bind must have at least ",
                     x@along, " dimensions for this binding operation"))

    dims <- get_dims_to_bind(x@seeds, x@along)
    if (is.character(dims))
        return(wmsg2(dims))
    TRUE
}

setValidity2("DelayedAbind", .validate_DelayedAbind)

new_DelayedAbind <- function(seeds, along)
{
    new2("DelayedAbind", seeds=seeds, along=along)
}

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
    Nindex <- vector(mode="list", length=length(index))
    Nindex[[x@along]] <- get_rev_index(part_idx)
    subset_by_Nindex(ans, Nindex)
}

setMethod("extract_array", "DelayedAbind", .extract_array_from_DelayedAbind)


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

