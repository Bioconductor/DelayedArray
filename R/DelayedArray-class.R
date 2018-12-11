### =========================================================================
### DelayedArray objects
### -------------------------------------------------------------------------


### The "root" node of the tree of DelayedOp objects. Represents a no-op.
setClass("DelayedArray", contains="DelayedUnaryIsoOp")

### Extending DataTable gives us a few things for free (head(), tail(),
### etc...). Note that even though DelayedMatrix already extends Array (via
### DelayedArray, DelayedUnaryOp, and DelayedOp) we need to make DelayedMatrix
### a *direct* child of Array and to place Array *before* DataTable in the
### 'contains' field below. This ensures that method dispatch will pick the
### method for Array in case a generic has methods defined for Array and
### DataTable (e.g. as.data.frame()). Furthermore, for some obscure reason,
### it seems that we also need to place all the classes that are in the
### inheritance path between DelayedArray and Array in the 'contains' field
### otherwise we get the following error when trying to instantiate a
### DelayedMatrix object with new("DelayedMatrix"):
###
###     Error: C stack usage  7971652 is too close to the limit
###
setClass("DelayedMatrix",
    contains=c("DelayedArray",
               "DelayedUnaryIsoOp", "DelayedUnaryOp", "DelayedOp",
               "Array",
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
### Check the internals of a DelayedArray object
###
### Internal representation of DelayedArray objects has changed in
### DelayedArray 0.5.11 (went from 5 slots to 3), then again in DelayedArray
### 0.5.24 (went from 3 slots to only 1):
###   - Internals v0: 5 slots (DelayedArray < 0.5.11)
###   - Internals v1: 3 slots (DelayedArray >= 0.5.11 and < 0.5.24)
###   - Internals v2: 1 slot  (DelayedArray >= 0.5.24)
###
### The helpers below detect the internals version used by an object. These
### helpers are used in the validity and "show" methods for DelayedArray
### objects.
###

.get_DelayedArray_internals_version <- function(object)
{
    if (.hasSlot(object, "metaindex") && .hasSlot(object, "is_transposed"))
        return(0L)
    if (.hasSlot(object, "index") && .hasSlot(object, "delayed_ops"))
        return(1L)
    return(2L)
}

.get_pkgversion_from_internals_version <- function(internals_version)
{
    switch(sprintf("v%d", internals_version),
           v0="< 0.5.11",
           v1=">= 0.5.11 and < 0.5.24",
           v2="current", # i.e. >= 0.5.24
           stop("invalid internals version"))
}

.get_DelayedArray_pkgversion <- function(object)
{
    internals_version <- .get_DelayedArray_internals_version(object)
    .get_pkgversion_from_internals_version(internals_version)
}

.msg_for_old_internals <- function(object, pkgversion)
{
    paste0(class(object), " object uses internal representation ",
           "from DelayedArray\n  ", pkgversion, " and cannot be displayed ",
           "or used. Please update it with:\n\n",
           "      object <- updateObject(object, verbose=TRUE)\n\n",
           "  and re-serialize it.")
}

.check_DelayedArray_internals <- function(object)
{
    pkgversion <- .get_DelayedArray_pkgversion(object)
    if (pkgversion != "current")
        stop(.msg_for_old_internals(object, pkgversion))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_DelayedArray <- function(x)
{
    pkgversion <- .get_DelayedArray_pkgversion(x)
    if (pkgversion != "current")
        return(paste0("\n  ", .msg_for_old_internals(x, pkgversion)))
#    seed_dim <- dim(x@seed)
#    seed_ndim <- length(seed_dim)
#    ## In the context of validObject(), 'class(x)' is always "DelayedArray"
#    ## and not the real class of 'x', which seems to be a bug in validObject().
#    ## This prevents us from doing the check below.
#    if (seed_ndim == 2L && !is(x, matrixClass(x)))
#        return(wmsg2("'x' has 2 dimensions but is not a ",
#                     matrixClass(x), " derivative"))
    TRUE
}

setValidity2("DelayedArray", .validate_DelayedArray)

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
### Seed contract
###
### We overwrite the methods for DelayedUnaryOp objects only to detect
### DelayedArray objects with old internals.
###

setMethod("dim", "DelayedArray",
    function(x)
    {
        .check_DelayedArray_internals(x)
        callNextMethod()
    }
)

setMethod("dimnames", "DelayedArray",
    function(x)
    {
        .check_DelayedArray_internals(x)
        callNextMethod()
    }
)

setMethod("extract_array", "DelayedArray",
    function(x, index)
    {
        .check_DelayedArray_internals(x)
        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Low-level constructor. Not intended to be used directly by the end user.
new_DelayedArray <- function(seed=new("array"), Class="DelayedArray")
{
    stopifnot(isSingleString(Class), extends(Class, "DelayedArray"))
    seed_ndim <- length(dim(seed))
    if (extends(Class, "DelayedMatrix")) {
        if (seed_ndim != 2L)
            stop(wmsg("the supplied seed must have exactly 2 dimensions ",
                      "when the specified class (", Class, ") extends ",
                      "DelayedMatrix"))
    } else {
        if (seed_ndim == 2L)
            Class <- matrixClass(new(Class))
    }
    new2(Class, seed=seed)
}

setGeneric("DelayedArray", function(seed) standardGeneric("DelayedArray"))

setMethod("DelayedArray", "ANY", function(seed) new_DelayedArray(seed))

setMethod("DelayedArray", "DelayedArray", function(seed) seed)

setMethod("DelayedArray", "DelayedOp",
    function(seed)
    {
        if (getOption("DelayedArray.simplify", default=TRUE)) {
            seed <- simplify(seed, incremental=TRUE)
            if (!is(seed, "DelayedOp"))
                return(DelayedArray(seed))
        }
        new_DelayedArray(seed)
    }
)


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
### simplify()
###

setMethod("simplify", "DelayedArray",
    function(x, incremental=FALSE)
        DelayedArray(simplify(x@seed, incremental=incremental))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### netSubsetAndAperm()
###

setMethod("netSubsetAndAperm", "DelayedArray",
    function(x, as.DelayedOp=FALSE)
    {
        x <- x@seed
        callGeneric()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### updateObject() must be able to update from internals v0 to v2 and from
### internals v1 to v2.
###

### Reflect internals v1.
setClass("DelayedArray1",
    contains="DelayedArray",
    representation(
        index="list",       # List (possibly named) of subscripts as
                            # positive integer vectors, one vector per
                            # seed dimension. *Missing* list elements
                            # are allowed and represented by NULLs.

        delayed_ops="list"  # List of delayed operations. See below
                            # for the details.
    ),
    prototype(
        seed=new("array"),
        index=list(NULL)
    )
)

.from_internals_v1_to_v2 <- function(object1)
{
    seed <- object1@seed
    seed_dimnames <- dimnames(seed)

    ## Translate 'index' slot as DelayedOp objects (1 DelayedSubset and
    ## 1 DelayedDimnames) and stash them inside 'seed'.

    index <- lapply(unname(object1@index), unname)
    op <- new2("DelayedSubset", seed=seed, index=index)
    if (!is_noop(op))
        seed <- op

    object1_dimnames <- lapply(seq_along(object1@index),
        function(along) {
            i <- object1@index[[along]]
            if (is.null(i)) seed_dimnames[[along]] else names(i)
        })
    if (all(S4Vectors:::sapply_isNULL(object1_dimnames)))
        object1_dimnames <- NULL

    op <- new_DelayedDimnames(seed, object1_dimnames)
    if (!is_noop(op))
        seed <- op

    ## Translate 'delayed_ops' slot as DelayedUnaryIsoOpWithArgs objects
    ## and stash them inside 'seed'.

    for (delayed_op in object1@delayed_ops) {
        OP <- delayed_op[[1L]]
        Largs <- delayed_op[[2L]]
        Rargs <- delayed_op[[3L]]
        Lalong <- Ralong <- NA
        recycle_along_last_dim <- delayed_op[[4L]]
        if (recycle_along_last_dim) {
            if (length(Largs) == 1L) Lalong <- 1L else Ralong <- 1L
        }
        seed <- new_DelayedUnaryIsoOpWithArgs(seed, OP, Largs, Rargs,
                                              Lalong, Ralong)
    }

    DelayedArray(seed)
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

.from_internals_v0_to_v2 <- function(object0)
{
    delayed_ops <- .update_delayed_ops(object0@delayed_ops)
    object1 <- new2("DelayedArray1", seed=object0@seed,
                                     index=object0@index,
                                     delayed_ops=delayed_ops,
                                     check=FALSE)
    object2 <- .from_internals_v1_to_v2(object1)

    if (identical(object0@metaindex, seq_along(object0@index)) &&
        identical(object0@is_transposed, FALSE))
        return(object2)

    if (any(vapply(delayed_ops,
                   function(delayed_op) delayed_op[[4L]],
                   logical(1))))
        stop(wmsg("object is too complex, sorry"))

    perm <- object0@metaindex
    if (object0@is_transposed)
        perm <- rev(perm)
    aperm(object2, perm)
}

.updateObject_DelayedArray <- function(object, ..., verbose=FALSE)
{
    object@seed <- updateObject(object@seed, verbose=verbose)

    internals_version <- .get_DelayedArray_internals_version(object)
    pkgversion <- .get_pkgversion_from_internals_version(internals_version)

    if (pkgversion == "current") {
        if (verbose)
            message("[updateObject] Internal representation of ",
                    class(object), " object is current.\n",
                    "[updateObject] Nothing to update.")
        return(object)
    }

    if (verbose)
        message("[updateObject] ", class(object), " object uses ",
                "internal representation from\n",
                "[updateObject] DelayedArray ", pkgversion, ". ",
                "Updating it ...")

    if (internals_version == 1L)
        return(.from_internals_v1_to_v2(object))

    if (internals_version == 0L)
        return(.from_internals_v0_to_v2(object))

    object
}

setMethod("updateObject", "DelayedArray", .updateObject_DelayedArray)


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
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedSubset(x@seed, Nindex)
    DelayedArray(op)
}

stash_DelayedAperm <- function(x, perm)
{
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedAperm(x@seed, perm)
    DelayedArray(op)
}

stash_DelayedUnaryIsoOpStack <- function(x, OP)
{
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedUnaryIsoOpStack(x@seed, OPS=list(OP), check.op=TRUE)
    DelayedArray(op)
}

stash_DelayedUnaryIsoOpWithArgs <- function(x, OP,
                                            Largs=list(), Rargs=list(),
                                            Lalong=NA, Ralong=NA)
{
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedUnaryIsoOpWithArgs(x@seed, OP=OP,
                                        Largs=Largs, Rargs=Rargs,
                                        Lalong=Lalong, Ralong=Ralong,
                                        check.op=TRUE)
    DelayedArray(op)
}

### stash_DelayedUnaryIsoOp() was introduced in BioC 3.7 and renamed
### stash_DelayedUnaryIsoOpWithArgs() in BioC 3.8. The alias below is for
### backward compatibility with objects that were serialized with BioC 3.7.
### For example, the 'BS.cancer.ex.fit' dataset (a BSseq instance from the
### bsseqData package and used in the bsseq_analysis.Rmd vignette from the
### bsseq package) contains the "plogis" method for DelayedArray objects in
### its "trans" slot, and this method contains a reference to
### stash_DelayedUnaryIsoOp. This is the sign that this BSseq instance was
### serialized with BioC 3.7. Without the alias below, the bsseq_analysis.Rmd
### vignette chokes on this line:
###
###   BS.cancer.ex.tstat <- BSmooth.tstat(BS.cancer.ex.fit, ...)
###
### with the following error:
###
###   Error: processing vignette 'bsseq_analysis.Rmd' failed with diagnostics:
###   could not find function "stash_DelayedUnaryIsoOp"
###
### Note that 'BS.cancer.ex.fit' contains 3 DelayedMatrix objects in
### its "assays" slot.
### TODO: Remove this alias at some point (e.g. in BioC 3.10).
stash_DelayedUnaryIsoOp <- stash_DelayedUnaryIsoOpWithArgs

stash_DelayedSubassign <- function(x, Nindex, value)
{
    stopifnot(is(x, "DelayedArray"))
    if (is(value, "DelayedArray"))
        value <- value@seed
    op <- new_DelayedSubassign(x@seed, Nindex, value)
    DelayedArray(op)
}

stash_DelayedDimnames <- function(x, dimnames)
{
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedDimnames(x@seed, dimnames)
    DelayedArray(op)
}

stash_DelayedAbind <- function(x, objects, along)
{
    stopifnot(is(x, "DelayedArray"))
    objects <- S4Vectors:::prepare_objects_to_bind(x, objects)
    seeds <- lapply(c(list(x), objects), slot, "seed")
    op <- new_DelayedAbind(seeds, along)
    DelayedArray(op)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Extend base::aperm() by allowing dropping and/or adding ineffective
### dimensions. See aperm2.R
###

### S3/S4 combo for aperm.DelayedArray
aperm.DelayedArray <- function(a, perm, ...) stash_DelayedAperm(a, perm, ...)
setMethod("aperm", "DelayedArray", aperm.DelayedArray)


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
        "be used to drop and/or add \"ineffective dimensions\" (i.e. ",
        "dimensions equal to 1) to the object."
    )
    can_map <- function() {
        if (length(idx1) != length(idx2))
            return(FALSE)
        if (length(idx1) == 0L)
            return(TRUE)
        all(new_dim[idx1] == old_dim[idx2])
    }
    if (!can_map())
        stop(wmsg(cannot_map_msg))

    new2old <- rep.int(NA_integer_, length(new_dim))
    new2old[idx1] <- idx2
    new2old
}

.set_DelayedArray_dim <- function(x, value)
{
    x_dim <- dim(x)
    value <- .normalize_dim_replacement_value(value, x_dim)
    new2old <- .map_new_to_old_dim(value, x_dim, class(x))
    ans <- aperm(x, new2old)
    stopifnot(identical(dim(ans), value))  # sanity check
    ans
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
### Show
###

setMethod("show", "DelayedArray",
    function(object)
    {
        .check_DelayedArray_internals(object)
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
    function(x, ..., recursive=FALSE)
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

