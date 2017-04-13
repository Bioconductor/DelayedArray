### =========================================================================
### DelayedArray objects
### -------------------------------------------------------------------------


setClass("DelayedArray",
    representation(
        seed="ANY",              # An array-like object expected to satisfy
                                 # the "seed contract" i.e. to support dim(),
                                 # dimnames(), and subset_seed_as_array().

        index="list",            # List (possibly named) of subscripts as
                                 # positive integer vectors, one vector per
                                 # seed dimension. *Missing* list elements
                                 # are allowed and represented by objects of
                                 # class "name".

        metaindex="integer",     # Index into the "index" slot specifying the
                                 # seed dimensions to keep.

        delayed_ops="list",      # List of delayed operations. See below
                                 # for the details.

        is_transposed="logical"  # Is the object considered to be transposed
                                 # with respect to the seed?
    ),
    prototype(
        seed=new("array"),
        index=list(integer(0)),
        metaindex=1L,
        is_transposed=FALSE
    )
)

### Extending DataTable gives us a few things for free (head(), tail(),
### etc...)
setClass("DelayedMatrix",
    contains=c("DelayedArray", "DataTable"),
    prototype=prototype(
        seed=new("matrix"),
        index=list(integer(0), integer(0)),
        metaindex=1:2
    )
)

### Overwrite unsafe automatic coercion method that return invalid objects (it
### doesn't validate them).
setAs("DelayedArray", "DelayedMatrix",
    function(from) new("DelayedMatrix", from)
)

### For internal use only.
setGeneric("matrixClass", function(x) standardGeneric("matrixClass"))

setMethod("matrixClass", "DelayedArray", function(x) "DelayedMatrix")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_DelayedArray <- function(x)
{
    x_ndim <- length(dim(x@seed))
    ## 'seed' slot.
    if (x_ndim == 0L)
        return(wmsg("'x@seed' must have dimensions"))
    ## 'index' slot.
    if (length(x@index) != x_ndim)
        return(wmsg("'x@index' must have one element per dimension ",
                    "in 'x@seed'"))
    ## 'metaindex' slot.
    if (length(x@metaindex) == 0L)
        return(wmsg("'x@metaindex' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@metaindex, 1L, x_ndim))
        return(wmsg("all values in 'x@metaindex' must be >= 1 ",
                    "and <= 'length(x@index)'"))
    if (!isStrictlySorted(x@metaindex))
        return(wmsg("'x@metaindex' must be strictly sorted"))
    if (!all(lengths(x@index)[-x@metaindex] == 1L))
        return(wmsg("all the dropped dimensions in 'x' must be equal to 1"))
    ## 'is_transposed' slot.
    if (!isTRUEorFALSE(x@is_transposed))
        return(wmsg("'x@is_transposed' must be TRUE or FALSE"))
    TRUE
}

setValidity2("DelayedArray", .validate_DelayedArray)

### TODO: Move this to S4Vectors and make it the validity method for DataTable
### object.
.validate_DelayedMatrix <- function(x)
{
    if (length(dim(x)) != 2L)
        return(wmsg("'x' must have exactly 2 dimensions"))
    TRUE
}

setValidity2("DelayedMatrix", .validate_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### NOT exported but used in HDF5Array!
new_DelayedArray <- function(seed=new("array"), Class="DelayedArray")
{
    seed <- remove_pristine_DelayedArray_wrapping(seed)
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    if (seed_ndim == 2L)
        Class <- matrixClass(new(Class))
    index <- rep.int(alist(foo=), seed_ndim)
    seed_dimnames <- dimnames(seed)
    ## If 'seed' has dimnames, then we put them on 'index'.
    if (!is.null(seed_dimnames)) {
        expand_idx <- which(!S4Vectors:::sapply_isNULL(seed_dimnames))
        index[expand_idx] <- mapply(
            function(d, names) setNames(seq_len(d), names),
            seed_dim[expand_idx], seed_dimnames[expand_idx],
            SIMPLIFY=FALSE, USE.NAMES=FALSE
        )
    }
    new2(Class, seed=seed, index=index, metaindex=seq_along(index))
}

DelayedArray <- function(seed) new_DelayedArray(seed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pristine objects
###
### A pristine DelayedArray object is an object that does not carry any
### delayed operation on it. In other words, it's in sync with (i.e. reflects
### the content of) its seed.
###

### NOT exported but used in HDF5Array!
### Note that false negatives happen when 'x' carries delayed operations that
### do nothing, but that's ok.
is_pristine <- function(x)
{
    ## 'x' should not carry any delayed operation on it, that is, all the
    ## DelayedArray slots must be in their original state.
    x2 <- new_DelayedArray(x@seed)
    class(x) <- class(x2) <- "DelayedArray"
    identical(x, x2)
}

### Remove the DelayedArray wrapping (or nested wrappings) from around the
### seed if the wrappings are pristine.
remove_pristine_DelayedArray_wrapping <- function(x)
{
    if (!(is(x, "DelayedArray") && is_pristine(x)))
        return(x)
    remove_pristine_DelayedArray_wrapping(x@seed)
}

### When a pristine DelayedArray derived object (i.e. an HDF5Array object) is
### about to be touched, we first need to downgrade it to a DelayedArray or
### DelayedMatrix *instance*.
downgrade_to_DelayedArray_or_DelayedMatrix <- function(x)
{
    if (is(x, "DelayedMatrix"))
        return(as(x, "DelayedMatrix"))
    as(x, "DelayedArray")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim()
###

### dim() getter.

.get_DelayedArray_dim_before_transpose <- function(x)
{
    get_subscripts_lengths(x@index, dim(x@seed))[x@metaindex]
}
.get_DelayedArray_dim <- function(x)
{
    ans <- .get_DelayedArray_dim_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}

setMethod("dim", "DelayedArray", .get_DelayedArray_dim)

### Even though prod() always returns a double, it seems that the length()
### primitive function automatically turns this double into an integer if
### it's <= .Machine$integer.max
setMethod("length", "DelayedArray", function(x) prod(dim(x)))

setMethod("isEmpty", "DelayedArray", function(x) any(dim(x) == 0L))

### dim() setter.

.normalize_dim_replacement_value <- function(value, x_dim)
{
    if (is.null(value))
        stop(wmsg("you can't do that, sorry"))
    if (!is.numeric(value))
        stop(wmsg("the supplied dim vector must be numeric"))
    if (length(value) == 0L)
        stop(wmsg("the supplied dim vector cannot be empty"))
    if (!is.integer(value))
        value <- as.integer(value)
    if (S4Vectors:::anyMissingOrOutside(value, 0L))
        stop(wmsg("the supplied dim vector cannot contain negative ",
                  "or NA values"))
    if (length(value) > length(x_dim))
        stop(wmsg("too many dimensions supplied"))
    prod1 <- prod(value)
    prod2 <- prod(x_dim)
    if (prod1 != prod2)
        stop(wmsg("the supplied dims [product ", prod1, "] do not match ",
                  "the length of object [", prod2, "]"))
    unname(value)
}

.map_new_to_old_dim <- function(new_dim, old_dim, x_class)
{
    idx1 <- which(new_dim != 1L)
    idx2 <- which(old_dim != 1L)

    cannot_map_msg <- wmsg(
        "Cannot map the supplied dim vector to the current dimensions of ",
        "the object. On a ", x_class, " object, the dim() setter can only ",
        "be used to drop some of the ineffective dimensions (the dimensions ",
        "equal to 1 are the ineffective dimensions)."
    )

    can_map <- function() {
        if (length(idx1) != length(idx2))
            return(FALSE)
        if (length(idx1) == 0L)
            return(TRUE)
        if (!all(new_dim[idx1] == old_dim[idx2]))
            return(FALSE)
        tmp <- idx2 - idx1
        tmp[[1L]] >= 0L && isSorted(tmp)
    }
    if (!can_map())
        stop(cannot_map_msg)

    new2old <- seq_along(new_dim) +
        rep.int(c(0L, idx2 - idx1), diff(c(1L, idx1, length(new_dim) + 1L)))

    if (new2old[[length(new2old)]] > length(old_dim))
        stop(cannot_map_msg)

    new2old
}

.set_DelayedArray_dim <- function(x, value)
{
    x_dim <- dim(x)
    value <- .normalize_dim_replacement_value(value, x_dim)
    new2old <- .map_new_to_old_dim(value, x_dim, class(x))
    stopifnot(identical(value, x_dim[new2old]))  # sanity check
    if (x@is_transposed) {
        x_metaindex <- rev(rev(x@metaindex)[new2old])
    } else {
        x_metaindex <- x@metaindex[new2old]
    }
    if (!identical(x@metaindex, x_metaindex)) {
        x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
        x@metaindex <- x_metaindex
    }
    if (length(dim(x)) == 2L)
        x <- as(x, matrixClass(x))
    x
}

setReplaceMethod("dim", "DelayedArray", .set_DelayedArray_dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

setMethod("drop", "DelayedArray",
    function(x)
    {
        x_dim <- dim(x)
        dim(x) <- x_dim[x_dim != 1L]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dimnames()
###

### dimnames() getter.

.get_DelayedArray_dimnames_before_transpose <- function(x)
{
    ans <- lapply(x@metaindex, function(N) names(x@index[[N]]))
    if (all(S4Vectors:::sapply_isNULL(ans)))
        return(NULL)
    ans
}
.get_DelayedArray_dimnames <- function(x)
{
    ans <- .get_DelayedArray_dimnames_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}

setMethod("dimnames", "DelayedArray", .get_DelayedArray_dimnames)

### dimnames() setter.

.normalize_dimnames_replacement_value <- function(value, ndim)
{
    if (is.null(value))
        return(vector("list", length=ndim))
    if (!is.list(value))
        stop("the supplied dimnames must be a list")
    if (length(value) > ndim)
        stop(wmsg("the supplied dimnames is longer ",
                  "than the number of dimensions"))
    if (length(value) <- ndim)
        length(value) <- ndim
    value
}

.set_DelayedArray_dimnames <- function(x, value)
{
    value <- .normalize_dimnames_replacement_value(value, length(x@metaindex))
    if (x@is_transposed)
        value <- rev(value)

    ## We quickly identify a no-op situation. While doing so, we are careful to
    ## not trigger a copy of the "index" slot (which can be big). The goal is
    ## to make a no-op like 'dimnames(x) <- dimnames(x)' as fast as possible.
    touched_midx <- which(mapply(
        function(N, names) !identical(names(x@index[[N]]), names),
        x@metaindex, value,
        USE.NAMES=FALSE
    ))
    if (length(touched_midx) == 0L)
        return(x)  # no-op

    x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
    touched_idx <- x@metaindex[touched_midx]
    x_seed_dim <- dim(x@seed)
    x@index[touched_idx] <- mapply(
        function(N, names) {
            i <- x@index[[N]]
            if (missing(i))
                i <- seq_len(x_seed_dim[[N]])  # expand 'i'
            setNames(i, names)
        },
        touched_idx, value[touched_midx],
        SIMPLIFY=FALSE, USE.NAMES=FALSE
    )
    x
}

setReplaceMethod("dimnames", "DelayedArray", .set_DelayedArray_dimnames)

### names() getter & setter.

.get_DelayedArray_names <- function(x)
{
    if (length(dim(x)) != 1L)
        return(NULL)
    dimnames(x)[[1L]]
}

setMethod("names", "DelayedArray", .get_DelayedArray_names)

.set_DelayedArray_names <- function(x, value)
{
    if (length(dim(x)) != 1L) {
        if (!is.null(value))
            stop("setting the names of a ", class(x), " object with more ",
                 "than 1 dimension is not supported")
        return(x)
    }
    dimnames(x)[[1L]] <- value
    x
}

setReplaceMethod("names", "DelayedArray", .set_DelayedArray_names)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [
###

### 'subscripts' must be a "multidimensional subscript" i.e. a list with one
### subscript per dimension in 'x'. Missing subscripts are represented by
### list elements of class "name".
.subset_DelayedArray_by_subscripts <- function(x, subscripts)
{
    x_index <- x@index
    x_ndim <- length(x@metaindex)
    x_seed_dim <- dim(x@seed)
    x_delayed_ops <- x@delayed_ops
    for (n in seq_along(subscripts)) {
        subscript <- subscripts[[n]]
        if (missing(subscript))
            next
        n0 <- if (x@is_transposed) x_ndim - n + 1L else n
        N <- x@metaindex[[n0]]
        i <- x_index[[N]]
        if (missing(i))
            i <- seq_len(x_seed_dim[[N]])  # expend 'i'
        x_index[[N]] <- extractROWS(i, subscript)
        if (n0 == 1L)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, subscript,
                                                      FALSE)
        if (n0 == x_ndim)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, subscript,
                                                      TRUE)
    }
    if (!identical(x@index, x_index)) {
        x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
        x@index <- x_index
        if (!identical(x@delayed_ops, x_delayed_ops))
            x@delayed_ops <- x_delayed_ops
    }
    x
}

### WARNING: This function uses the ugly "substitute() + eval()" hack! One
### problem with this hack is that it seems to be causing the following bug:
###
###     library(DelayedArray)
###     toto <- function(a, z) a[ , , z]
###     a <- array(1:60, c(3, 5, 4))
###     A <- DelayedArray(a)
###     toto(a, 4L)
###     toto(A, 4L)
###     # Error in eval(expr, envir, enclos) : object 'z' not found
###
### This is probably due to the somewhat arbitrary use of
### 'envir=parent.frame(2L)' in the call to eval().
### FIXME: The above bug needs to be fixed!
.extract_DelayedArray_subset <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")

    ## Check the dimensionality of the user call i.e whether the function was
    ## called with 1D-style, or 2D-style, or 3D-style etc... subsetting.
    ndim <- nargs() - 1L
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (!missing(drop))
        ndim <- ndim - 1L
    if (ndim == 1L && missing(i))
        ndim <- 0L
    if (ndim != 0L && ndim != x_ndim) {
        if (ndim == 1L)
            stop(wmsg("1D-style subsetting is not supported on DelayedArray ",
                      "objects (yet)"))
        stop("incorrect number of dimensions")
    }

    ## Prepare the "multidimensional subscript".
    subscripts <- rep.int(alist(foo=), x_ndim)
    if (!missing(i))
        subscripts[[1L]] <- i
    if (!missing(j))
        subscripts[[2L]] <- j
    ## Ugly "substitute() + eval()" hack! See above.
    dots <- substitute(...())  # list of non-evaluated args
    for (n2 in seq_along(dots)) {
        k <- dots[[n2]]
        if (!missing(k))
            subscripts[[2L + n2]] <- eval(k, envir=parent.frame(2L))
    }

    ## Perform the subsetting.
    .subset_DelayedArray_by_subscripts(x, subscripts)
}

setMethod("[", "DelayedArray", .extract_DelayedArray_subset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()
###

### Return the slice as a list.
.extract_data_frame_slice <- function(x, index)
{
    slice <- subset_by_subscripts(x, index)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}
.extract_DataFrame_slice <- function(x, index)
{
    slice <- subset_by_subscripts(x, index)
    slice <- as.data.frame(slice)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}

### Return a list with one list element per column in data frame 'x'.
### All the list elements have length 0.
.extract_data_frame_slice0 <- function(x)
{
    slice0 <- x[0L, , drop=FALSE]
    ## Turn into a list and replace factors with character vectors.
    lapply(slice0, as.vector)
}
.extract_DataFrame_slice0 <- function(x)
{
    slice0 <- x[0L, , drop=FALSE]
    slice0 <- as.data.frame(slice0)
    if (ncol(slice0) != ncol(x))
        stop(wmsg("DataFrame object 'x' can be used as the seed of ",
                  "a DelayedArray object only if as.data.frame(x) ",
                  "preserves the number of columns"))
    ## Turn into a list and replace factors with character vectors.
    lapply(slice0, as.vector)
}

### Equivalent to 'typeof(as.matrix(x))' but with an almost-zero
### memory footprint (it avoids the cost of turning 'x' into a matrix).
.get_data_frame_type <- function(x)
{
    if (ncol(x) == 0L)
        return("logical")
    slice0 <- .extract_data_frame_slice0(x)
    typeof(unlist(slice0, use.names=FALSE))
}

### Equivalent to 'typeof(as.matrix(as.data.frame(x)))' but with an
### almost-zero memory footprint (it avoids the cost of turning 'x' first
### into a data frame then into a matrix).
.get_DataFrame_type <- function(x)
{
    if (ncol(x) == 0L)
        return("logical")
    slice0 <- .extract_DataFrame_slice0(x)
    typeof(unlist(slice0, use.names=FALSE))
}

### 'index' is expected to be an unnamed list of subscripts as positive integer
### vectors, one vector per seed dimension. *Missing* list elements are allowed
### and represented by objects of class "name".
### The "subset_seed_as_array" methods don't need to support anything else.
### They must return an ordinary array. No need to propagate the dimnames.
setGeneric("subset_seed_as_array", signature="seed",
    function(seed, index) standardGeneric("subset_seed_as_array")
)

setMethod("subset_seed_as_array", "ANY",
    function(seed, index)
    {
        slice <- subset_by_subscripts(seed, index)
        as.array(slice)
    }
)

setMethod("subset_seed_as_array", "array",
    function(seed, index)
        subset_by_subscripts(seed, index)
)

### Equivalent to
###
###     subset_by_subscripts(as.matrix(x), index)
###
### but avoids the cost of turning the full data frame 'x' into a matrix so
### memory footprint stays small when 'index' is small.
setMethod("subset_seed_as_array", "data.frame",
    function(seed, index)
    {
        #ans_type <- .get_data_frame_type(seed)
        slice0 <- .extract_data_frame_slice0(seed)
        slice <- .extract_data_frame_slice(seed, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_subscripts_lengths(index, dim(seed)))
    }
)

### Equivalent to
###
###     subset_by_subscripts(as.matrix(as.data.frame(x)), index)
###
### but avoids the cost of turning the full DataFrame 'x' first into a data
### frame then into a matrix so memory footprint stays small when 'index' is
### small.
setMethod("subset_seed_as_array", "DataFrame",
    function(seed, index)
    {
        #ans_type <- .get_DataFrame_type(seed)
        slice0 <- .extract_DataFrame_slice0(seed)
        slice <- .extract_DataFrame_slice(seed, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_subscripts_lengths(index, dim(seed)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Management of delayed operations
###
### The 'delayed_ops' slot represents the list of delayed operations op1, op2,
### etc... Each delayed operation is itself represented by a list of length 4:
###   1) The name of the function to call (e.g. "+" or "log").
###   2) The list of "left arguments" i.e. the list of arguments to place
###      before the array in the function call.
###   3) The list of "right arguments" i.e. the list of arguments to place
###      after the array in the function call.
###   4) A single logical. Indicates the dimension along which the (left or
###      right) argument of the function call needs to be recycled when the
###      operation is actually executed (done by .execute_delayed_ops() which
###      is called by as.array()). FALSE: along the 1st dim; TRUE: along
###      the last dim; NA: no recycling. Recycling is only supported for
###      function calls with 2 arguments (i.e. the array and the recycled
###      argument) at the moment.
###      
### Each operation must return an array of the same dimensions as the original
### array.
###

register_delayed_op <- function(x, FUN, Largs=list(), Rargs=list(),
                                        recycle_along_last_dim=NA)
{
    if (isTRUEorFALSE(recycle_along_last_dim)) {
        nLargs <- length(Largs)
        nRargs <- length(Rargs)
        ## Recycling is only supported for function calls with 2 arguments
        ## (i.e. the array and the recycled argument) at the moment.
        stopifnot(nLargs + nRargs == 1L)
        partially_recycled_arg <- if (nLargs == 1L) Largs[[1L]] else Rargs[[1L]]
        stopifnot(length(partially_recycled_arg) == nrow(x))
    }
    delayed_op <- list(FUN, Largs, Rargs, recycle_along_last_dim)
    x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
    x@delayed_ops <- c(x@delayed_ops, list(delayed_op))
    x
}

.subset_delayed_op_args <- function(delayed_op, i, subset_along_last_dim)
{
    recycle_along_last_dim <- delayed_op[[4L]]
    if (is.na(recycle_along_last_dim)
     || recycle_along_last_dim != subset_along_last_dim)
        return(delayed_op)
    Largs <- delayed_op[[2L]]
    Rargs <- delayed_op[[3L]]
    nLargs <- length(Largs)
    nRargs <- length(Rargs)
    stopifnot(nLargs + nRargs == 1L)
    if (nLargs == 1L) {
        new_arg <- extractROWS(Largs[[1L]], i)
        delayed_op[[2L]] <- list(new_arg)
    } else {
        new_arg <- extractROWS(Rargs[[1L]], i)
        delayed_op[[3L]] <- list(new_arg)
    }
    if (length(new_arg) == 1L)
        delayed_op[[4L]] <- NA
    delayed_op
}

.subset_delayed_ops_args <- function(delayed_ops, i, subset_along_last_dim)
    lapply(delayed_ops, .subset_delayed_op_args, i, subset_along_last_dim)

### 'a' is the ordinary array returned by the "combining" operator.
.execute_delayed_ops <- function(a, delayed_ops)
{
    a_dim <- dim(a)
    first_dim <- a_dim[[1L]]
    last_dim <- a_dim[[length(a_dim)]]
    a_len <- length(a)
    if (a_len == 0L) {
        p1 <- p2 <- 0L
    } else {
        p1 <- a_len / first_dim
        p2 <- a_len / last_dim
    }

    recycle_arg <- function(partially_recycled_arg, recycle_along_last_dim) {
        if (recycle_along_last_dim) {
            stopifnot(length(partially_recycled_arg) == last_dim)
            rep(partially_recycled_arg, each=p2)
        } else {
            stopifnot(length(partially_recycled_arg) == first_dim)
            rep.int(partially_recycled_arg, p1)
        }
    }

    prepare_call_args <- function(a, delayed_op) {
        Largs <- delayed_op[[2L]]
        Rargs <- delayed_op[[3L]]
        recycle_along_last_dim <- delayed_op[[4L]]
        if (isTRUEorFALSE(recycle_along_last_dim)) {
            nLargs <- length(Largs)
            nRargs <- length(Rargs)
            stopifnot(nLargs + nRargs == 1L)
            if (nLargs == 1L) {
                Largs <- list(recycle_arg(Largs[[1L]], recycle_along_last_dim))
            } else {
                Rargs <- list(recycle_arg(Rargs[[1L]], recycle_along_last_dim))
            }
        }
        c(Largs, list(a), Rargs)
    }

    for (delayed_op in delayed_ops) {
        FUN <- delayed_op[[1L]]
        call_args <- prepare_call_args(a, delayed_op)

        ## Perform the delayed operation.
        a <- do.call(FUN, call_args)

        ## Some vectorized operations on an ordinary array can drop the dim
        ## attribute (e.g. comparing a zero-col matrix with an atomic vector).
        a_new_dim <- dim(a)
        if (is.null(a_new_dim)) {
            ## Restore the dim attribute.
            dim(a) <- a_dim
        } else {
            ## Sanity check.
            stopifnot(identical(a_dim, a_new_dim))
        }
    }
    a
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transpose
###

### The actual transposition of the data is delayed i.e. it will be realized
### on the fly only when as.array() (or as.vector() or as.matrix()) is called
### on 'x'.
setMethod("t", "DelayedArray",
    function(x)
    {
        x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
        x@is_transposed <- !x@is_transposed
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### as.array()
###

### TODO: Not sure we need this. Using drop() should do it.
.reduce_array_dimensions <- function(x)
{
    x_dim <- dim(x)
    x_dimnames <- dimnames(x)
    effdim_idx <- which(x_dim != 1L)  # index of effective dimensions
    if (length(effdim_idx) >= 2L) {
        dim(x) <- x_dim[effdim_idx]
        dimnames(x) <- x_dimnames[effdim_idx]
    } else {
        dim(x) <- NULL
        if (length(effdim_idx) == 1L)
            names(x) <- x_dimnames[[effdim_idx]]
    }
    x
}

### Realize the object i.e. execute all the delayed operations and turn the
### object back into an ordinary array.
.from_DelayedArray_to_array <- function(x, drop=FALSE)
{
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    a <- subset_seed_as_array(x@seed, unname(x@index))
    dim(a) <- .get_DelayedArray_dim_before_transpose(x)
    ans <- .execute_delayed_ops(a, x@delayed_ops)
    dimnames(ans) <- .get_DelayedArray_dimnames_before_transpose(x)
    if (drop)
        ans <- .reduce_array_dimensions(ans)
    ## Base R doesn't support transposition of an array of arbitrary dimension
    ## (generalized transposition) so the call to t() below will fail if 'ans'
    ## has more than 2 dimensions. If we want as.array() to work on a
    ## transposed DelayedArray object of arbitrary dimension, we need to
    ## implement our own generalized transposition of an ordinary array.
    if (x@is_transposed) {
        if (length(dim(ans)) > 2L)
            stop("can't do as.array() on this object, sorry")
        ans <- t(ans)
    }
    ans
}

### S3/S4 combo for as.array.DelayedArray
as.array.DelayedArray <- function(x, ...)
    .from_DelayedArray_to_array(x, ...)
setMethod("as.array", "DelayedArray", .from_DelayedArray_to_array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercions based on as.array()
###

### S3/S4 combo for as.vector.DelayedArray
as.vector.DelayedArray <- function(x, mode="any")
{
    ans <- as.array(x, drop=TRUE)
    as.vector(ans, mode=mode)
}
setMethod("as.vector", "DelayedArray", as.vector.DelayedArray)

slicing_tip <- c(
    "Consider reducing its number of effective dimensions by slicing it ",
    "first (e.g. x[8, 30, , 2, ]). Make sure that all the indices used for ",
    "the slicing have length 1 except at most 2 of them which can be of ",
    "arbitrary length or missing."
)

.from_DelayedArray_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (sum(x_dim != 1L) > 2L)
        stop(wmsg(class(x), " object with more than 2 effective dimensions ",
                  "cannot be coerced to a matrix. ", slicing_tip))
    ans <- as.array(x, drop=TRUE)
    if (length(x_dim) == 2L) {
        dim(ans) <- x_dim
        dimnames(ans) <- dimnames(x)
    } else {
        as.matrix(ans)
    }
    ans
}

### S3/S4 combo for as.matrix.DelayedArray
as.matrix.DelayedArray <- function(x, ...) .from_DelayedArray_to_matrix(x, ...)
setMethod("as.matrix", "DelayedArray", .from_DelayedArray_to_matrix)

### S3/S4 combo for as.data.frame.DelayedArray
as.data.frame.DelayedArray <- function(x, row.names=NULL, optional=FALSE, ...)
    as.data.frame(as.array(x, drop=TRUE),
                  row.names=row.names, optional=optional, ...)
setMethod("as.data.frame", "DelayedArray", as.data.frame.DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type()
###
### For internal use only.
###

setGeneric("type", function(x) standardGeneric("type"))

setMethod("type", "array", function(x) typeof(x))

### If 'x' is a DelayedArray object, 'type(x)' must always return the same
### as 'typeof(as.array(x))'.
setMethod("type", "DelayedArray",
    function(x)
    {
        subscripts <- as.list(integer(length(dim(x))))
        ## x0 <- x[0, ..., 0]
        x0 <- .subset_DelayedArray_by_subscripts(x, subscripts)
        typeof(as.array(x0, drop=TRUE))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [[
###

.get_DelayedArray_element <- function(x, i)
{
    i <- normalizeDoubleBracketSubscript(i, x)
    subscripts <- as.integer(arrayInd(i, dim(x)))
    as.vector(.subset_DelayedArray_by_subscripts(x, subscripts))
}

### Only support linear subscripting at the moment.
### TODO: Support multidimensional subscripting e.g. x[[5, 15, 2]] or
### x[["E", 15, "b"]].
setMethod("[[", "DelayedArray",
    function(x, i, j, ...)
    {
        dots <- list(...)
        if (length(dots) > 0L)
            dots <- dots[names(dots) != "exact"]
        if (!missing(j) || length(dots) > 0L)
            stop("incorrect number of subscripts")
        .get_DelayedArray_element(x, i)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "DelayedArray", show_compact_array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining and splitting
###
### Combining arrays with c() is NOT an endomorphism!
###

### 'objects' must be a list of array-like objects that support as.vector().
combine_array_objects <- function(objects)
{
    if (!is.list(objects))
        stop("'objects' must be a list")

    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    if (length(objects) == 0L)
        return(NULL)

    unlist(lapply(objects, as.vector), recursive=FALSE, use.names=FALSE)
}

setMethod("c", "DelayedArray",
    function (x, ..., recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for DelayedArray objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
        } else {
            objects <- list(x, ...)
        }
        combine_array_objects(objects)
    }
)

setMethod("splitAsList", "DelayedArray",
    function(x, f, drop=FALSE, ...)
        splitAsList(as.vector(x), f, drop=drop, ...)
)

### S3/S4 combo for split.DelayedArray
split.DelayedArray <- function(x, f, drop=FALSE, ...)
    splitAsList(x, f, drop=drop, ...)
setMethod("split", c("DelayedArray", "ANY"), split.DelayedArray)

