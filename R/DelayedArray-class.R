### =========================================================================
### DelayedArray objects
### -------------------------------------------------------------------------


### The "root" node of the tree of DelayedOp objects. Represents a no-op.
setClass("DelayedArray", contains="DelayedUnaryOp")

### Extending DataTable gives us a few things for free (head(), tail(),
### etc...). Note that even though DelayedMatrix already extends Array (via
### DelayedArray, DelayedUnaryOp, and DelayedOp) we need to make DelayedMatrix
### a *direct* child of Array and to place Array *before* DataTable in the
### 'contains' field below. This ensures that method dispatch will pick the
### method for Array in case a generic has methods defined for Array and
### DataTable (e.g. as.data.frame()). Furthermore, it seems that we also need
### to place all the classes that are in the inheritance path between
### DelayedArray and Array in the 'contains' field otherwise we get the
### following error when trying to instantiate a DelayedMatrix object with
### new("DelayedMatrix"):
###
###     Error: C stack usage  7971652 is too close to the limit
###
setClass("DelayedMatrix",
    contains=c("DelayedArray", "DelayedUnaryOp", "DelayedOp", "Array",
               "DataTable"),
    prototype=prototype(
        seed=new("matrix")
    )
)

### Automatic coercion method from DelayedArray to DelayedMatrix silently
### returns a broken object (unfortunately these dummy automatic coercion
### methods don't bother to validate the object they return). So we overwrite
### it.
setAs("DelayedArray", "DelayedMatrix",
    function(from) new("DelayedMatrix", from)
)

### The user should not be able to degrade a DelayedMatrix object to a
### DelayedArray object so 'as(x, "DelayedArray", strict=TRUE)' should
### fail or be a no-op when 'x' is a DelayedMatrix object. Making this
### coercion a no-op seems to be the easiest (and safest) way to go.
setAs("DelayedMatrix", "DelayedArray", function(from) from)  # no-op

### For internal use only.
setGeneric("matrixClass", function(x) standardGeneric("matrixClass"))

setMethod("matrixClass", "DelayedArray", function(x) "DelayedMatrix")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#.validate_DelayedArray <- function(x)
#{
#    seed_dim <- dim(x@seed)
#    seed_ndim <- length(seed_dim)
#    ## In the context of validObject(), 'class(x)' is always "DelayedArray"
#    ## and not the real class of 'x', which seems to be a bug in validObject().
#    ## This prevents us from doing the check below.
#    if (seed_ndim == 2L && !is(x, matrixClass(x)))
#        return(wmsg2("'x' has 2 dimensions but is not a ",
#                     matrixClass(x), " derivative"))
#    TRUE
#}
#
#setValidity2("DelayedArray", .validate_DelayedArray)

### TODO: Move this to S4Vectors and make it the validity method for DataTable
### object.
.validate_DelayedMatrix <- function(x)
{
    if (length(dim(x)) != 2L)
        return(wmsg2("'x' must have exactly 2 dimensions"))
    TRUE
}

setValidity2("DelayedMatrix", .validate_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### NOT exported but used in HDF5Array!
new_DelayedArray <- function(seed=new("array"), Class="DelayedArray")
{
    seed_ndim <- length(dim(seed))
    if (seed_ndim == 2L)
        Class <- matrixClass(new(Class))
    new2(Class, seed=seed)
}

setGeneric("DelayedArray", function(seed) standardGeneric("DelayedArray"))

setMethod("DelayedArray", "ANY",
    function(seed)
    {
        if (getOption("DelayedArray.simplify", default=TRUE))
            seed <- simplify(seed)
        new_DelayedArray(seed)
    }
)

setMethod("DelayedArray", "DelayedArray", function(seed) seed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### Internal representation of DelayedArray objects has changed in
### DelayedArray 0.5.11 (Bioc 3.7).
###

.get_DelayedArray_version <- function(object)
{
    if (.hasSlot(object, "metaindex") && .hasSlot(object, "is_transposed"))
        return("< 0.5.11")
    if (.hasSlot(object, "index") && .hasSlot(object, "delayed_ops"))
        return(">= 0.5.11 and < 0.5.24")
    "current"  # i.e. >= 0.5.24
}

.update_delayed_ops <- function(delayed_ops)
{
    for (i in seq_along(delayed_ops)) {
        delayed_op <- delayed_ops[[i]]
        recycle_along_last_dim <- delayed_op[[4L]]
        if (is.na(recycle_along_last_dim)) {
            recycle_along_first_dim <- FALSE
        } else if (!recycle_along_last_dim) {
            recycle_along_first_dim <- TRUE
        } else {
            stop(wmsg("object is too complex, sorry"))
        }
        delayed_op[[4L]] <- recycle_along_first_dim
        delayed_ops[[i]] <- delayed_op
    }
    delayed_ops
}

.updateObject_DelayedArray <- function(object, ..., verbose=FALSE)
{
    object@seed <- updateObject(object@seed, verbose=verbose)

    version <- .get_DelayedArray_version(object)

    if (version == "current") {
        if (verbose)
            message("[updateObject] Internal representation of ",
                    class(object), " object is current.\n",
                    "[updateObject] Nothing to update.")
        return(object)
    }

    if (verbose)
        message("[updateObject] ", class(object), " object uses ",
                "internal representation from\n",
                "[updateObject] DelayedArray ", version, ". Updating it ...")

    if (version == "< 0.5.11") {
        delayed_ops <- .update_delayed_ops(object@delayed_ops)
        ## THIS DOES NOT WORK ANYMORE! (no more 'index' or 'delayed_ops' slot
        ## in DelayedArray >= 0.5.24)
        ans <- new2(class(object), seed=object@seed,
                                   index=object@index,
                                   delayed_ops=delayed_ops,
                                   check=FALSE)
        if (identical(object@metaindex, seq_along(object@index)) &&
            identical(object@is_transposed, FALSE))
            return(ans)
        if (any(vapply(delayed_ops,
                       function(delayed_op) delayed_op[[4L]],
                       logical(1))))
            stop(wmsg("object is too complex, sorry"))
        perm <- object@metaindex
        if (object@is_transposed)
            perm <- rev(perm)
        return(aperm(ans, perm))
    }

    if (version == ">= 0.5.11 and < 0.5.24") {
        seed <- object@seed
        seed_dimnames <- dimnames(seed)

        ## Translate 'index' slot as DelayedOp objects (1 DelayedSubset and
        ## 1 DelayedDimnames) and stash them inside 'seed'.

        index <- lapply(unname(object@index), unname)
        op <- new2("DelayedSubset", seed=seed, index=index)
        if (!isNoOp(op))
            seed <- op

        object_dimnames <- lapply(seq_along(object@index),
            function(along) {
                i <- object@index[[along]]
                if (is.null(i)) seed_dimnames[[along]] else names(i)
            })
        if (all(S4Vectors:::sapply_isNULL(object_dimnames)))
            object_dimnames <- NULL

        op <- new_DelayedDimnames(seed, object_dimnames)
        if (!isNoOp(op))
            seed <- op

        ## Translate 'delayed_ops' slot as DelayedUnaryIsoOp objects and
        ## stash them inside 'seed'.

        for (delayed_op in object@delayed_ops) {
            OP <- delayed_op[[1L]]
            Largs <- delayed_op[[2L]]
            Rargs <- delayed_op[[3L]]
            Lalong <- Ralong <- NA
            recycle_along_last_dim <- delayed_op[[4L]]
            if (recycle_along_last_dim) {
                if (length(Largs) == 1L) Lalong <- 1L else Ralong <- 1L
            }
            seed <- new_DelayedUnaryIsoOp(seed, OP, Largs, Rargs,
                                                Lalong, Ralong)
        }

        return(DelayedArray(seed))
    }

    object
}

setMethod("updateObject", "DelayedArray", .updateObject_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simplify()
###

setMethod("simplify", "DelayedArray",
    function(x) DelayedArray(simplify(x@seed))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarizeMappingToSeed()
###

setMethod("summarizeMappingToSeed", "DelayedArray",
    function(x)
    {
        x <- x@seed
        callGeneric()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Stash delayed ops in a DelayedArray object
###
### Each "stashing" utility below:
###   1) creates a new DelayedOp object representing a particular delayed
###      operation,
###   2) "stashes" the new DelayedOp object inside the input DelayedArray
###      object 'x' by inserting it in the tree of DelayedOp objects (stored
###      in 'x@seed') as the new root node.
###
### No "stashing" utility for nodes of type DelayedNaryIsoOp for now.
### See DelayedOp-class.R for the list of all node types.
###

stash_DelayedSubset <- function(x, Nindex)
{
    op <- new_DelayedSubset(x@seed, Nindex)
    DelayedArray(op)
}

stash_DelayedAperm <- function(x, perm)
{
    op <- new_DelayedAperm(x@seed, perm)
    DelayedArray(op)
}

stash_DelayedUnaryIsoOp <- function(x, OP, Largs=list(), Rargs=list(),
                                       Lalong=NA, Ralong=NA)
{
    op <- new_DelayedUnaryIsoOp(x@seed, OP=OP, Largs=Largs, Rargs=Rargs,
                                        Lalong=Lalong, Ralong=Ralong)
    DelayedArray(op)
}

stash_DelayedDimnames <- function(x, dimnames)
{
    op <- new_DelayedDimnames(x@seed, dimnames)
    DelayedArray(op)
}

stash_DelayedAbind <- function(..., along)
{
    seeds <- lapply(unname(list(...)), slot, "seed")
    op <- new_DelayedAbind(seeds, along)
    DelayedArray(op)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pristine objects
###
### A pristine DelayedArray object is an object that does not carry any
### delayed operation.
###

### Note that false negatives happen when 'x' carries delayed operations that
### do nothing, but that's ok.
is_pristine <- function(x) { !is(x@seed, "DelayedOp") }


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nseed()
###
### Return the number of leaves in the tree of DelayedOp objects contained
### in a DelayedArray object.
###

setGeneric("nseed", function(x) standardGeneric("nseed"))

setMethod("nseed", "ANY",
    function(x)
    {
        if (is(x, "DelayedUnaryOp"))
            return(nseed(x@seed))
        if (is(x, "DelayedNaryOp"))
            x <- x@seeds
        if (is.list(x)) {
            ans <- sum(vapply(x, nseed, integer(1), USE.NAMES=FALSE))
            return(ans)
        }
        return(1L)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seed() getter/setter
###

### If the tree of DelayedOp objects contained in the DelayedArray object
### has a single leaf (i.e. if the tree is linear), then seed() returns it.
### Otherwise, it raises an error.
setGeneric("seed", function(x) standardGeneric("seed"))

setMethod("seed", "DelayedOp",
    function(x)
    {
        if (is(x, "DelayedNaryOp")) {
            ## Tree is not linear.
            stop(wmsg("seed() ", IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS))
        }
        x1 <- x@seed
        if (!is(x1, "DelayedOp"))
            return(x1)  ## found the leaf seed
        x <- x1
        callGeneric()  # recursive call
    }
)

### If the tree of DelayedOp objects contained in the DelayedArray object
### has a single leaf (i.e. if the tree is linear), then the seed() setter
### replaces it. Otherwise, it raises an error.
setGeneric("seed<-", signature="x",
    function(x, ..., value) standardGeneric("seed<-")
)

.normalize_seed_replacement_value <- function(value, x_seed)
{
    if (!is(value, class(x_seed)))
        stop(wmsg("supplied seed must be a ", class(x_seed), " object"))
    if (!identical(dim(value), dim(x_seed)))
        stop(wmsg("supplied seed must have the same dimensions ",
                  "as current seed"))
    if (!identical(dimnames(value), dimnames(x_seed)))
        stop(wmsg("supplied seed must have the same dimnames ",
                  "as current seed"))
    value
}

setReplaceMethod("seed", "DelayedOp",
    function(x, value)
    {
        if (is(x, "DelayedNaryOp")) {
            ## Tree is not linear.
            stop(wmsg("the seed() setter ", IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS))
        }
        x1 <- x@seed
        if (!is(x1, "DelayedOp")) {
            ## Replace the leaf seed.
            x@seed <- .normalize_seed_replacement_value(value, x1)
            return(x)
        }
        seed(x@seed) <- value  # recursive call
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### The path of a DelayedArray object is the path of its leaf seed. So path()
### will work on a DelayedArray object only if it works on its leaf seed.
### For example it will work if its leaf seed is an on-disk object (e.g. an
### HDF5ArraySeed object) but not if it's an in-memory object (e.g. an
### ordinary array or RleArraySeed object).
setMethod("path", "DelayedArray",
    function(object, ...) path(seed(object), ...)
)

### The path() setter will work on a DelayedArray object only if it works on
### its leaf seed. For example it will work if its leaf seed is an on-disk
### object (e.g. an HDF5ArraySeed object) but not if it's an in-memory object
### (e.g. an ordinary array or RleArraySeed object).
setReplaceMethod("path", "DelayedArray",
    function(object, ..., value)
    {
        path(seed(object), ...) <- value
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Unlike base::aperm(), the method for DelayedArray objects supports
### dropping dimensions. Note that only "ineffective" dimensions can be
### dropped (i.e. dimensions equal to 1, so dropping them preserves the
### length). This feature is used by the "drop" method below.
###

setGeneric("aperm", signature="a")

.aperm.DelayedArray <- function(a, perm)
{
    if (missing(perm))
        perm <- rev(seq_along(dim(a)))
    stash_DelayedAperm(a, perm)
}

### S3/S4 combo for aperm.DelayedArray
aperm.DelayedArray <- function(a, perm, ...) .aperm.DelayedArray(a, perm, ...)
setMethod("aperm", "DelayedArray", aperm.DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

setMethod("drop", "DelayedArray",
    function(x)
    {
        perm <- which(dim(x) != 1L)
        if (length(perm) == 0L)
            perm <- 1L
        aperm(x, perm)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() setter
###
### On a DelayedArray object, the dim() setter can only be used to drop some
### of the "ineffective" dimensions (i.e. dimensions equal to 1, so dropping
### them preserves the length of the object).
###

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

    cannot_map_msg <- c(
        "Cannot map the supplied dim vector to the current dimensions of ",
        "the object. On a ", x_class, " object, the dim() setter can only ",
        "be used to drop some of the \"ineffective\" dimensions (i.e. ",
        "dimensions equal to 1, so dropping them preserves the length of ",
        "the object)."
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
        stop(wmsg(cannot_map_msg))

    new2old <- seq_along(new_dim) +
        rep.int(c(0L, idx2 - idx1), diff(c(1L, idx1, length(new_dim) + 1L)))

    if (new2old[[length(new2old)]] > length(old_dim))
        stop(wmsg(cannot_map_msg))

    new2old
}

.set_DelayedArray_dim <- function(x, value)
{
    x_dim <- dim(x)
    value <- .normalize_dim_replacement_value(value, x_dim)
    new2old <- .map_new_to_old_dim(value, x_dim, class(x))
    stopifnot(identical(value, x_dim[new2old]))  # sanity check
    aperm(x, new2old)
}

setReplaceMethod("dim", "DelayedArray", .set_DelayedArray_dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dimnames() setter
###

setReplaceMethod("dimnames", "DelayedArray",
    function(x, value) stash_DelayedDimnames(x, value)
)

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
            stop(wmsg("setting the names of a ", class(x), " object ",
                      "with more than 1 dimension is not supported"))
        return(x)
    }
    dimnames(x)[[1L]] <- value
    x
}

setReplaceMethod("names", "DelayedArray", .set_DelayedArray_names)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [
###

### Linear single bracket subsetting e.g. x[5].
### NOT delayed (i.e. return an atomic vector).
.subset_DelayedArray_by_1Dindex <- function(x, i)
{
    x_dim <- dim(x)
    i_dim <- dim(i)
    if (identical(x_dim, i_dim) && type(i) == "logical") {
        i <- which(i)  # calls .DelayedArray_block_which() defined
                       # in DelayedArray-utils.R
    } else {
        i <- normalizeSingleBracketSubscript2(i, length(x))
    }
    i_len <- length(i)
    if (i_len == 0L) {
        ## x0 <- x[integer(0), ..., integer(0)]
        index <- rep.int(list(integer(0)), length(x_dim))
        x0 <- extract_array(x, index)
        return(as.vector(x0))
    }
    if (i_len == 1L)
        return(.get_DelayedArray_element(x, i))

    ## We want to walk only on the blocks that we actually need to visit so we
    ## don't use block_APPLY() or family because they walk on all the blocks.

    grid <- defaultGrid(x)
    nblock <- length(grid)

    breakpoints <- cumsum(lengths(grid))
    part_idx <- get_part_index(i, breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    block_idx <- which(lengths(split_part_idx) != 0L)  # blocks to visit
    res <- lapply(block_idx, function(b) {
            if (get_verbose_block_processing())
                message("Visiting block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            block <- extract_block(x, grid[[b]])
            if (!is.array(block))
                block <- as.array(block)
            block_ans <- block[split_part_idx[[b]]]
            if (get_verbose_block_processing())
                message("OK")
            block_ans
    })
    unlist(res, use.names=FALSE)[get_rev_index(part_idx)]
}

.subset_DelayedArray <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(x)  # no-op
    x_ndim <- length(dim(x))
    if (nsubscript == 1L && (x_ndim != 1L || drop)) {
        ## Linear single bracket subsetting e.g. x[5].
        ## If 'x' is mono-dimensional and 'drop' is FALSE, we switch
        ## to "multi-dimensional single bracket subsetting" which is
        ## delayed.
        return(.subset_DelayedArray_by_1Dindex(x, Nindex[[1L]]))
    }
    if (nsubscript != x_ndim)
        stop("incorrect number of subscripts")
    ## Multi-dimensional single bracket subsetting
    ##
    ##     x[i_1, i_2, ..., i_n]
    ##
    ## Return an object of the same class as 'x' (endomorphism).
    ans <- stash_DelayedSubset(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "DelayedArray", .subset_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [<-
###

.filling_error_msg <- c(
    "filling a DelayedArray object 'x' with a value (i.e. 'x[] <- value') ",
    "is supported only when 'value' is an atomic vector and 'length(value)' ",
    "is a divisor of 'nrow(x)'"
)

.subassign_error_msg <- c(
    "subassignment to a DelayedArray object 'x' (i.e. 'x[i] <- value') is ",
    "supported only when the subscript 'i' is a logical DelayedArray object ",
    "with the same dimensions as 'x' and when 'value' is a scalar (i.e. an ",
    "atomic vector of length 1)"
)

.fill_DelayedArray_with_value <- function(x, value)
{
    if (!(is.vector(value) && is.atomic(value)))
        stop(wmsg(.filling_error_msg))
    value_len <- length(value)
    if (value_len == 1L)
        return(stash_DelayedUnaryIsoOp(x, `[<-`, Rargs=list(value=value)))
    x_len <- length(x)
    if (value_len > x_len)
        stop(wmsg("'value' is longer than 'x'"))
    x_nrow <- nrow(x)
    if (x_nrow != 0L) {
        if (value_len == 0L || x_nrow %% value_len != 0L)
            stop(wmsg(.filling_error_msg))
        value <- rep(value, length.out=x_nrow)
    }
    stash_DelayedUnaryIsoOp(x, `[<-`, Rargs=list(value=value), Ralong=1L)
}

.subassign_DelayedArray <- function(x, i, j, ..., value)
{
    if (missing(x))
        stop("'x' is missing")
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(.fill_DelayedArray_with_value(x, value))
    if (nsubscript != 1L)
        stop(wmsg(.subassign_error_msg))
    i <- Nindex[[1L]]
    if (!(is(i, "DelayedArray") &&
          identical(dim(x), dim(i)) &&
          type(i) == "logical"))
        stop(wmsg(.subassign_error_msg))
    if (!(is.vector(value) && is.atomic(value) && length(value) == 1L))
        stop(wmsg(.subassign_error_msg))
    DelayedArray(new_DelayedNaryIsoOp(x@seed, i@seed,
                                      OP=`[<-`, Rargs=list(value=value)))
}

setReplaceMethod("[", "DelayedArray", .subassign_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to sparse matrix (requires the Matrix package)
###

.from_DelayedMatrix_to_dgCMatrix <- function(from)
{
    idx <- which(from != 0L)
    array_ind <- arrayInd(idx, dim(from))
    i <- array_ind[ , 1L]
    j <- array_ind[ , 2L]
    x <- from[idx]
    Matrix::sparseMatrix(i, j, x=x, dims=dim(from), dimnames=dimnames(from))
}

setAs("DelayedMatrix", "dgCMatrix", .from_DelayedMatrix_to_dgCMatrix)
setAs("DelayedMatrix", "sparseMatrix", .from_DelayedMatrix_to_dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [[
###
### TODO: Implement a "getArrayElement" method for Array objects (in
### extract_array.R) that is based on extract_array() e.g it just needs
### to do something like:
###
###     as.vector(extract_array(x, as.list(subscripts)))
###
### This will make multi-dimensional and linear [[ work on DelayedArray
### objects. Then remove the method below and update DelayedArray-class.Rd
###

.get_DelayedArray_element <- function(x, i)
{
    i <- normalizeDoubleBracketSubscript(i, x)
    index <- as.list(arrayInd(i, dim(x)))
    as.vector(extract_array(x, index))
}

### Only support linear subsetting at the moment.
### TODO: Support multi-dimensional subsetting e.g. x[[5, 15, 2]] or
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

setMethod("show", "DelayedArray",
    function(object)
    {
        version <- .get_DelayedArray_version(object)
        if (version != "current")
            stop(class(object), " object uses internal representation from ",
                 "DelayedArray ", version, "\n  and cannot be displayed or ",
                 "used. Please update it with:\n",
                 "    object <- updateObject(object, verbose=TRUE)")
        show_compact_array(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining and splitting
###
### TODO: These methods are based on as.vector() so could be made methods
### for Array objects (and moved to extract_array.R).
###

### Note that combining arrays with c() is NOT an endomorphism!
setMethod("c", "DelayedArray",
    function (x, ..., recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop(wmsg("\"c\" method for DelayedArray objects ",
                      "does not support the 'recursive' argument"))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Binding
###
### We only support binding DelayedArray objects along the rows or the cols
### at the moment. No binding along an arbitrary dimension yet! (i.e. no
### "abind" method yet)
###

### arbind() and acbind()

.DelayedArray_arbind <- function(...) stash_DelayedAbind(..., along=1L)
.DelayedArray_acbind <- function(...) stash_DelayedAbind(..., along=2L)

setMethod("arbind", "DelayedArray", .DelayedArray_arbind)
setMethod("acbind", "DelayedArray", .DelayedArray_acbind)

### rbind() and cbind()

setMethod("rbind", "DelayedMatrix", .DelayedArray_arbind)
setMethod("cbind", "DelayedMatrix", .DelayedArray_acbind)

.as_DelayedMatrix_objects <- function(objects)
{
    lapply(objects,
        function(object) {
            if (length(dim(object)) != 2L)
                stop(wmsg("cbind() and rbind() are not supported on ",
                          "DelayedArray objects that don't have exactly ",
                          "2 dimensions. Please use acbind() or arbind() ",
                          "instead."))
            as(object, "DelayedMatrix")
        })
}

.DelayedArray_rbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("rbind", objects)
}

.DelayedArray_cbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("cbind", objects)
}

setMethod("rbind", "DelayedArray", .DelayedArray_rbind)
setMethod("cbind", "DelayedArray", .DelayedArray_cbind)

