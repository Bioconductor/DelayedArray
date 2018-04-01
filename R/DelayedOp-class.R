### =========================================================================
### DelayedOp objects
### -------------------------------------------------------------------------
###
### In a DelayedArray object the delayed operations are stored as a tree of
### DelayedOp objects. Each node in the tree is represented by a DelayedOp
### object. 5 types of nodes are currently supported. Each type is a concrete
### DelayedOp subclass:
###
###   Node type      outdegree  operation
###   ---------------------------------------------------------------------
###   DelayedSubset          1  Multi-dimensional single bracket subsetting
###   DelayedUnaryIsoOp      1  Unary op that preserves the geometry
###   DelayedAperm           1  Extended aperm() (can drop dimensions)
###   DelayedVariadicIsoOp   N  N-ary op that preserves the geometry
###   DelayedAbind           N  abind()
###
### All the nodes are array-like objects that must satisfy the "seed contract"
### i.e. they must support dim(), dimnames(), and extract_array().
### Unary nodes (i.e. nodes with an outdegree of 1) must also support the
### seed() and path() getters and setters.
###

### This virtual class and its 5 concrete subclasses are for internal use
### only and are not exported.
setClass("DelayedOp", contains="Array", representation("VIRTUAL"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Default seed() and path() getters/setters
###

setGeneric("seed", function(x) standardGeneric("seed"))

setGeneric("seed<-", signature="x",
    function(x, ..., value) standardGeneric("seed<-")
)

.IS_NOT_SUPOORTED_ETC <- c(
    "is not supported on a DelayedArray object with ",
    "multiple leaf seeds at the moment"
)

setMethod("seed", "DelayedOp",
    function(x)
    {
        if (.hasSlot(x, "seeds"))
            stop(wmsg("seed() ", .IS_NOT_SUPOORTED_ETC))
        x@seed
    }
)

normalize_seed_replacement_value <- function(value, x_seed)
{
    if (!is(value, class(x_seed)))
        stop(wmsg("supplied seed must be a ", class(x_seed), " object"))
    if (!identical(dim(value), dim(x_seed)))
        stop(wmsg("supplied seed must have the same dimensions ",
                  "as current seed"))
    value
}

setReplaceMethod("seed", "DelayedOp",
    function(x, value)
    {
        if (.hasSlot(x, "seeds"))
            stop(wmsg("the seed() setter ", .IS_NOT_SUPOORTED_ETC))
        x@seed <- normalize_seed_replacement_value(value, seed(x))
        x
    }
)

setMethod("path", "DelayedOp",
    function(object, ...)
    {
        if (.hasSlot(object, "seeds"))
            stop(wmsg("path() ", .IS_NOT_SUPOORTED_ETC))
        ## The path of a DelayedOp object is the path of its seed so path()
        ## will work only on a DelayedOp object that has a seed that supports
        ## path().
        ## For example it will work if the seed is an on-disk object (e.g. an
        ## HDF5ArraySeed object) but not if it's an in-memory object (e.g. an
        ## ordinary array or RleArraySeed object).
        path(seed(object), ...)
    }
)

setReplaceMethod("path", "DelayedOp",
    function(object, ..., value)
    {
        if (.hasSlot(object, "seeds"))
            stop(wmsg("the path() setter ", .IS_NOT_SUPOORTED_ETC))
        ## The path() setter sets the supplied path on the seed of the
        ## DelayedOp object so it will work out-of-the-box on any DelayedOp
        ## object that has a seed that supports the path() setter.
        ## For example it will work if the seed is an HDF5ArraySeed object.
        path(seed(object), ...) <- value
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedSubset objects
###

setClass("DelayedSubset",
    contains="DelayedOp",
    representation(
        seed="ANY",   # An array-like object expected to satisfy the "seed
                      # contract".

        index="list"  # List (possibly named) of subscripts as positive
                      # integer vectors, one per seed dimension. *Missing*
                      # list elements are allowed and represented by NULLs.
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

### Seed contract.

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
### DelayedUnaryIsoOp objects
###

setClass("DelayedUnaryIsoOp",
    contains="DelayedOp",
    representation(
        seed="ANY",     # An array-like object expected to satisfy the
                        # "seed contract".

        OP="function",  # The function to apply to the seed (e.g. `+` or
                        # log). It should act as an isomorphism i.e. always
                        # return an array parallel to the input array (i.e.
                        # same dimensions).

        Largs="list",   # Additional left arguments to OP.

        Rargs="list"    # Additional right arguments to OP.
    ),
    prototype(
        seed=new("array"),
        OP=identity
    )
)

new_DelayedUnaryIsoOp <- function(seed=new("array"),
                                  OP=identity, Largs=list(), Rargs=list())
{
    OP <- match.fun(OP)
    new2("DelayedUnaryIsoOp", seed=seed, OP=OP, Largs=Largs, Rargs=Rargs)
}

### Seed contract.

setMethod("dim", "DelayedUnaryIsoOp", function(x) dim(x@seed))

setMethod("dimnames", "DelayedUnaryIsoOp", function(x) dimnames(x@seed))

setMethod("extract_array", "DelayedUnaryIsoOp",
    function(x, index)
    {
        a <- extract_array(x@seed, index)
        do.call(x@OP, c(x@Largs, list(a), x@Rargs))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedAperm objects
###

setClass("DelayedAperm",
    contains="DelayedOp",
    representation(
        seed="ANY",                # An array-like object expected to satisfy
                                   # the "seed contract".

        dim_combination="integer"  # Index into dim(seed) specifying the seed
                                   # dimensions to keep.
    ),
    prototype(
        seed=new("array"),
        dim_combination=1L
    )
)

.validate_DelayedAperm <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)
    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))
    ## 'dim_combination' slot.
    if (length(x@dim_combination) == 0L)
        return(wmsg2("'x@dim_combination' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@dim_combination, 1L, seed_ndim))
        return(wmsg2("all values in 'x@dim_combination' must be >= 1 ",
                     "and <= 'seed_ndim'"))
    if (anyDuplicated(x@dim_combination))
        return(wmsg2("'x@dim_combination' cannot have duplicates"))
    if (!all(seed_dim[-x@dim_combination] == 1L))
        return(wmsg2("dimensions to drop from 'x' must be equal to 1"))
    TRUE
}

setValidity2("DelayedAperm", .validate_DelayedAperm)

new_DelayedAperm <- function(seed, dim_combination=NULL)
{
    if (is.null(dim_combination))
        dim_combination <- seq_along(dim(seed))
    new2("DelayedAperm", seed=seed, dim_combination=dim_combination)
}

### Seed contract.

.get_DelayedAperm_dim <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_dim[x@dim_combination]
}

setMethod("dim", "DelayedAperm", .get_DelayedAperm_dim)

.get_DelayedAperm_dimnames <- function(x)
{
    seed_dimnames <- dimnames(x@seed)
    seed_dimnames[x@dim_combination]  # return NULL if 'seed_dimnames' is NULL
}

setMethod("dimnames", "DelayedAperm", .get_DelayedAperm_dimnames)

.extract_array_from_DelayedAperm <- function(x, index)
{
    seed_dim <- dim(x@seed)
    seed_index <- rep.int(list(1L), length(seed_dim))
    seed_index[x@dim_combination] <- index
    a <- extract_array(x@seed, seed_index)
    dim(a) <- dim(a)[sort(x@dim_combination)]
    aperm(a, perm=rank(x@dim_combination))
}

setMethod("extract_array", "DelayedAperm",
    .extract_array_from_DelayedAperm
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DelayedVariadicIsoOp objects
###

setClass("DelayedVariadicIsoOp",
    contains="DelayedOp",
    representation(
        seeds="list",   # List of conformable array-like objects to combine.
                        # Each object is expected to satisfy the "seed
                        # contract".

        OP="function",  # The function to use to combine the seeds. It should
                        # act as an isomorphism i.e. always return an array
                        # parallel to the input arrays (i.e. same dimensions).

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

### Seed contract.

setMethod("dim", "DelayedVariadicIsoOp", function(x) dim(x@seeds[[1L]]))

setMethod("dimnames", "DelayedVariadicIsoOp",
    function(x) combine_dimnames(x@seeds)
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

setClass("DelayedAbind",
    contains="DelayedOp",
    representation(
        seeds="list",    # List of array-like objects to bind. Each object
                         # is expected to satisfy the "seed contract".

        along="integer"  # Single integer indicating the dimension along
                         # which to bind the seeds.
    ),
    prototype(
        seeds=list(new("array")),
        along=1L
    )
)

.validate_DelayedAbind <- function(x)
{
    if (length(x@seeds) == 0L)
        return(wmsg2("'x@seeds' cannot be empty"))
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

