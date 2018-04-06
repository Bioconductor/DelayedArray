### =========================================================================
### DelayedOp objects
### -------------------------------------------------------------------------
###
### In a DelayedArray object the delayed operations are stored as a tree of
### DelayedOp objects. Each node in the tree is represented by a DelayedOp
### object. 6 types of nodes are currently supported. Each type is a concrete
### DelayedOp subclass:
###
###   Node type      Outdegree  Operation
###   ---------------------------------------------------------------------
###   DelayedSubset          1  Multi-dimensional single bracket subsetting
###   DelayedDimnames        1  Set dimnames
###   DelayedUnaryIsoOp      1  Unary op that preserves the geometry
###   DelayedAperm           1  Extended aperm() (can drop dimensions)
###   DelayedVariadicIsoOp   N  N-ary op that preserves the geometry
###   DelayedAbind           N  abind()
###
### All the nodes are array-like objects that must comply with the "seed
### contract" i.e. they must support dim(), dimnames(), and extract_array().
###

### This virtual class and its 6 concrete subclasses are for internal use
### only and are not exported.
setClass("DelayedOp", contains="Array", representation("VIRTUAL"))

### NOT exported for now.
setGeneric("isNoOp", function(x) standardGeneric("isNoOp"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### showtree()
###
### A much more condensed version of str().
###

setGeneric("showtree", signature="x",
    function(x, show.node.dim=TRUE) standardGeneric("showtree")
)

### Avoid use of non-ASCII characters in R source code. There must be a much
### better way to do this.
.VBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x82)))
.TEE   <- rawToChar(as.raw(c(0xe2, 0x94, 0x9c)))
.ELBOW <- rawToChar(as.raw(c(0xe2, 0x94, 0x94)))
.HBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x80)))

### 'last.child' can be NA, TRUE, or FALSE. NA means 'x' is the root of the
### tree.
.show_tree <- function(x, indent="", last.child=NA, show.node.dim=TRUE)
{
    stopifnot(isSingleString(indent))
    stopifnot(is.logical(last.child), length(last.child) == 1L)

    if (!is.list(x)) {

        ## Display summary line.

        if (is.na(last.child)) {
            ## No prefix.
            prefix <- ""
        } else {
            ## 3-char prefix
            prefix <- paste0(if (last.child) .ELBOW else .TEE, .HBAR, " ")
        }

        is_leaf <- FALSE
        if (is(x, "DelayedOp")) {
            x_as1string <- summary(x)
        } else {
            x_as1string <- paste0(class(x), " object")
            if (!(.hasSlot(x, "seed") || .hasSlot(x, "seeds"))) {
                is_leaf <- TRUE
                x_as1string <- paste0("[seed] ", x_as1string)
            }
        }
        if (show.node.dim) {
            dim_in1string <- paste0(dim(x), collapse="x")
            x_as1string <- sprintf("%s %s: %s", dim_in1string, type(x),
                                                x_as1string)
        }
        cat(indent, prefix, x_as1string, "\n", sep="")
        if (is_leaf)
            return(invisible(NULL))
    }

    ## Display children.

    if (!is.na(last.child)) {
        ## Increase indent by 3 chars.
        indent <- paste0(indent, if (last.child) " " else .VBAR, "  ")
    }
    if (.hasSlot(x, "seed"))
        .show_tree(x@seed, indent, last.child=TRUE,
                   show.node.dim=show.node.dim)
    if (.hasSlot(x, "seeds"))
        x <- x@seeds
    if (is.list(x)) {
        nseed <- length(x)
        for (i in seq_len(nseed))
            .show_tree(x[[i]], indent, last.child=(i==nseed),
                       show.node.dim=show.node.dim)
    }
}

setMethod("showtree", "ANY",
    function(x, show.node.dim=TRUE)
    {
        if (!isTRUEorFALSE(show.node.dim))
            stop("'show.node.dim' must be TRUE or FALSE")
        .show_tree(x, show.node.dim=show.node.dim)
    }
)

setMethod("show", "DelayedOp", function(object) .show_tree(object))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedSubset objects
###
### Delayed "Multi-dimensional single bracket subsetting".
###

setClass("DelayedSubset",
    contains="DelayedOp",
    representation(
        seed="ANY",   # The input array-like object. Expected to comply with
                      # the "seed contract".

        index="list"  # List of subscripts as positive integer vectors, one
                      # per dimension in the input. *Missing* list elements
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

### 'Nindex' must be a "multidimensional subsetting Nindex" (see utils.R).
new_DelayedSubset <- function(seed=new("array"), Nindex=list(NULL))
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    stopifnot(is.list(Nindex), length(Nindex) == seed_ndim)

    ## Normalize 'Nindex' i.e. check and turn its non-NULL list elements into
    ## positive integer vectors.
    seed_dimnames <- dimnames(seed)
    index <- lapply(seq_len(seed_ndim),
                    function(along) {
                        subscript <- Nindex[[along]]
                        if (is.null(subscript))
                            return(NULL)
                        normalizeSingleBracketSubscript2(subscript,
                                                         seed_dim[[along]],
                                                         seed_dimnames[[along]])
                    })
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
    function(x) subset_dimnames(dimnames(x@seed), x@index)
)

.extract_array_from_DelayedSubset <- function(x, index)
{
    seed_dim <- dim(x@seed)
    stopifnot(is.list(index), length(index) == length(seed_dim))
    index2 <- lapply(seq_along(x@index),
                     function(along) {
                         i1 <- x@index[[along]]
                         i2 <- index[[along]]
                         if (is.null(i2))
                             return(i1)
                         if (is.null(i1))
                             return(i2)
                         i1[i2]
                     })
    extract_array(x@seed, index2)
}

setMethod("extract_array", "DelayedSubset", .extract_array_from_DelayedSubset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedDimnames objects
###
### Delayed "Set dimnames".
###

.INHERIT_FROM_SEED <- -1L

setClass("DelayedDimnames",
    contains="DelayedOp",
    representation(
        seed="ANY",      # The input array-like object. Expected to comply
                         # with the "seed contract".

        dimnames="list"  # List with one list element per dimension in
                         # the input. Each list element must be NULL,
                         # or a character vector, or special value
                         # .INHERIT_FROM_SEED
    ),
    prototype(
        seed=new("array"),
        dimnames=list(.INHERIT_FROM_SEED)
    )
)

.validate_DelayedDimnames <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)

    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))

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

setMethod("dim", "DelayedDimnames", function(x) dim(x@seed))

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

setMethod("extract_array", "DelayedDimnames",
    function(x, index) extract_array(x@seed, index)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedUnaryIsoOp objects
###
### Delayed "Unary op that preserves the geometry".
###

setClass("DelayedUnaryIsoOp",
    contains="DelayedOp",
    representation(
        seed="ANY",      # The input array-like object. Expected to comply
                         # with the "seed contract".

        OP="function",   # The function to apply to the input (e.g. `+` or
                         # log). It should act as an isomorphism i.e. always
                         # return an array-like object *parallel* to the input
                         # (i.e. with the same dimensions as the input).

        Largs="list",    # Left arguments to OP i.e. arguments to place
                         # before the input array in the function call.
        Rargs="list",    # Right arguments to OP i.e. arguments to place
                         # after the input array in the function call.

        Lidx="integer",  # Index of left arguments that are parallel to the
                         # rows of the input.
        Ridx="integer"   # Index of right arguments that are parallel to the
                         # rows of the input.
    ),
    prototype(
        seed=new("array"),
        OP=identity
    )
)

new_DelayedUnaryIsoOp <- function(seed=new("array"),
                                  OP=identity, Largs=list(), Rargs=list(),
                                  Lidx=integer(0), Ridx=integer(0))
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    OP <- match.fun(OP)

    seed_nrow <- seed_dim[[1L]]
    stopifnot(is.list(Largs), is.list(Rargs),
              is.integer(Lidx), is.integer(Ridx),
              all(elementNROWS(Largs[Lidx]) == seed_nrow),
              all(elementNROWS(Rargs[Ridx]) == seed_nrow))

    new2("DelayedUnaryIsoOp", seed=seed, OP=OP, Largs=Largs, Rargs=Rargs,
                              Lidx=Lidx, Ridx=Ridx)
}

### S3/S4 combo for summary.DelayedUnaryIsoOp

.DelayedUnaryIsoOp_summary <- function(object) "Unary op"

summary.DelayedUnaryIsoOp <-
    function(object, ...) .DelayedUnaryIsoOp_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOp", summary.DelayedUnaryIsoOp)

### Seed contract.

setMethod("dim", "DelayedUnaryIsoOp", function(x) dim(x@seed))

setMethod("dimnames", "DelayedUnaryIsoOp", function(x) dimnames(x@seed))

setMethod("extract_array", "DelayedUnaryIsoOp",
    function(x, index)
    {
        a <- extract_array(x@seed, index)
        Largs <- x@Largs
        Rargs <- x@Rargs
        i1 <- index[[1L]]
        if (!is.null(i1)) {
            Largs[x@Lidx] <- lapply(Largs[x@Lidx], extractROWS, i1)
            Rargs[x@Ridx] <- lapply(Rargs[x@Ridx], extractROWS, i1)
        }
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
### DelayedAperm objects
###
### Delayed "Extended aperm()" (can drop dimensions).
### Note that only "ineffective" dimensions can be dropped (i.e. dimensions
### equal to 1, so dropping them preserves the length).
###

setClass("DelayedAperm",
    contains="DelayedOp",
    representation(
        seed="ANY",     # The input array-like object. Expected to comply
                        # with the "seed contract".

        perm="integer"  # Index into dim(seed) specifying the *rearrangement*
                        # of the dimensions i.e. which dimensions of the input
                        # to keep and in which order.
    ),
    prototype(
        seed=new("array"),
        perm=1L
    )
)

### The seed is referred to as 'a' in the error messages. This makes them
### more meaningful for the end user in the context of calling aperm().
.validate_DelayedAperm <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)

    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))

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

.normarg_perm <- function(perm, a_dim)
{
    if (is.null(perm))
        return(seq_along(a_dim))
    if (!is.numeric(perm))
        stop(wmsg("'perm' must be an integer vector"))
    if (!is.integer(perm))
        perm <- as.integer(perm)
    perm
}

new_DelayedAperm <- function(seed, perm=NULL)
{
    perm <- .normarg_perm(perm, dim(seed))
    new2("DelayedAperm", seed=seed, perm=perm)
}

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
### DelayedVariadicIsoOp objects
###
### Delayed "N-ary op that preserves the geometry".
###

setClass("DelayedVariadicIsoOp",
    contains="DelayedOp",
    representation(
        seeds="list",   # The input array-like objects. Each object is
                        # expected to comply with the "seed contract".
                        # The objects must be "conformable" i.e. they must
                        # all have the same dimensions.

        OP="function",  # The function to use to combine the input objects.
                        # Should act as an isomorphism i.e. always return an
                        # array-like object *parallel* to the input objects
                        # (i.e. with the same dimensions).

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        seeds=list(new("array")),
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

.validate_DelayedVariadicIsoOp <- function(x)
{
    ## 'seeds' slot.
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))
    if (!.objects_are_conformable_arrays(x@seeds))
        return(wmsg2("'x@seeds' must be a list of conformable ",
                     "array-like objects"))
    TRUE
}

setValidity2("DelayedVariadicIsoOp", .validate_DelayedVariadicIsoOp)

new_DelayedVariadicIsoOp <- function(seed=new("array"), ...,
                                     OP=identity, Rargs=list())
{
    seeds <- unname(list(seed, ...))
    OP <- match.fun(OP)
    new2("DelayedVariadicIsoOp", seeds=seeds, OP=OP, Rargs=Rargs)
}

### S3/S4 combo for summary.DelayedVariadicIsoOp

.DelayedVariadicIsoOp_summary <- function(object) "N-ary op"

summary.DelayedVariadicIsoOp <-
    function(object, ...) .DelayedVariadicIsoOp_summary(object, ...)

setMethod("summary", "DelayedVariadicIsoOp", summary.DelayedVariadicIsoOp)

### Seed contract.

setMethod("dim", "DelayedVariadicIsoOp", function(x) dim(x@seeds[[1L]]))

setMethod("dimnames", "DelayedVariadicIsoOp",
    function(x) simplify_NULL_dimnames(combine_dimnames(x@seeds))
)

setMethod("extract_array", "DelayedVariadicIsoOp",
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
    contains="DelayedOp",
    representation(
        seeds="list",    # The input array-like objects. Each object is
                         # expected to comply with the "seed contract".

        along="integer"  # Single integer indicating the dimension along
                         # which to bind the input objects.
    ),
    prototype(
        seeds=list(new("array")),
        along=1L
    )
)

.validate_DelayedAbind <- function(x)
{
    ## 'seeds' slot.
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))

    ## 'along' slot.
    if (!(isSingleInteger(x@along) && x@along > 0L))
        return(wmsg2("'x@along' must be a single positive integer"))

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

### S3/S4 combo for summary.DelayedVariadicIsoOp

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
### SeedBinder classes were renamed DelayedAperm, DelayedVariadicIsoOp, and
### DelayedAbind, respectively.
### DelayedArray objects serialized with DelayedArray < 0.5.24 might contain
### instances of these old classes nested in their "seed" slot so we need to
### keep them around for now.
###

setClass("SeedDimPicker", contains="DelayedAperm")
setClass("ConformableSeedCombiner", contains="DelayedVariadicIsoOp")
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
        object <- new2("DelayedVariadicIsoOp", seeds=object@seeds,
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

