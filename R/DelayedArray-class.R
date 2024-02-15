### =========================================================================
### DelayedArray objects
### -------------------------------------------------------------------------


### The "root" node of the tree of DelayedOp objects. Represents a no-op.
setClass("DelayedArray", contains="DelayedUnaryIsoOp")

### Extending RectangularData gives us a few things for free (e.g. validity
### method for RectangularData objects, head(), tail(), etc...). Note
### that even though DelayedMatrix already extends Array (via DelayedArray,
### DelayedUnaryOp, and DelayedOp) we need to make DelayedMatrix a *direct*
### child of Array and to place Array *before* RectangularData in
### the 'contains' field below. This ensures that method dispatch will pick
### the method for Array in case a generic has methods defined for Array and
### RectangularData. Furthermore, for some obscure reason, it seems that we
### also need to place all the classes that are in the inheritance path
### between DelayedArray and Array in the 'contains' field otherwise we get
### the following error when trying to instantiate a DelayedMatrix object
### with new("DelayedMatrix"):
###
###     Error: C stack usage  7971652 is too close to the limit
###
setClass("DelayedMatrix",
    contains=c("DelayedArray",
               "DelayedUnaryIsoOp", "DelayedUnaryOp", "DelayedOp",
               "Array",
               "RectangularData"),
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
           "from DelayedArray\n  ", pkgversion, " and cannot be ",
           "displayed or used. Please update it with:\n\n",
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

### Note that, because DelayedArray extends DelayedUnaryOp, validation of
### a DelayedArray object is mostly taken care of by the validity method
### for DelayedUnaryOp objects, which is defined in DelayedOp-class.R.
.validate_DelayedArray <- function(x)
{
    pkgversion <- .get_DelayedArray_pkgversion(x)
    if (pkgversion != "current")
        return(.msg_for_old_internals(x, pkgversion))
#    seed_dim <- dim(x@seed)
#    seed_ndim <- length(seed_dim)
#    ## In the context of validObject(), 'class(x)' is always "DelayedArray"
#    ## and not the real class of 'x', which seems to be a bug in validObject().
#    ## This prevents us from doing the check below.
#    if (seed_ndim == 2L && !is(x, matrixClass(x)))
#        return(paste0("'x' has 2 dimensions but is not a ",
#                      matrixClass(x), " derivative"))
    TRUE
}

setValidity2("DelayedArray", .validate_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###
### We overwrite the methods for DelayedUnaryOp objects for the sole purpose
### of detecting DelayedArray objects with old internals.
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
    ## 'dimnames(seed)' can fail e.g. if 'seed' is or contains an
    ## HDF5ArraySeed object that points to a non-existing file, but we still
    ## want to be able to update object1.
    ## Our use case for this is ExperimentHub resource EH1656. This is a
    ## SummarizedExperiment object (added to ExperimentHub on 2017-10-06
    ## by the restfulSEData folks) where the assay is a very old DelayedMatrix
    ## instance (predates DelayedArray 0.4) that binds together 14 old
    ## HDF5ArraySeed instances that point to a non-existing file ('assays.h5').
    seed_dimnames <- try(dimnames(seed), silent=TRUE)
    if (inherits(seed_dimnames, "try-error"))
        seed_dimnames <- NULL

    ## Translate 'index' slot as DelayedOp objects (1 DelayedSubset and
    ## 1 DelayedSetDimnames) and stash them inside 'seed'.

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

    op <- new_DelayedSetDimnames(seed, object1_dimnames)
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
    object@seed <- updateObject(object@seed, ..., verbose=verbose)

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

stash_DelayedSetDimnames <- function(x, dimnames)
{
    stopifnot(is(x, "DelayedArray"))
    op <- new_DelayedSetDimnames(x@seed, dimnames)
    DelayedArray(op)
}

stash_DelayedAbind <- function(x, objects, along)
{
    stopifnot(is.list(objects))
    all_objects <- unname(S4Vectors:::delete_NULLs(c(list(x), objects)))
    ## All supplied objects must be DelayedArray objects or derivatives.
    seeds <- lapply(all_objects,
        function(object) {
            if (!is(object, "DelayedArray"))
                stop(wmsg("all the supplied objects must be DelayedArray ",
                          "objects or derivatives"))
            object@seed
        })
    op <- new_DelayedAbind(seeds, along)
    DelayedArray(op)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Like aperm2() in S4Arrays, extend base::aperm() by allowing dropping
### and/or adding ineffective dimensions.
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

### Return the mapping as an integer vector 'new2old' parallel to 'new_dim'.
### 'new2old' has the following properties:
### 1. It can contain NA's but it must have at least one non-NA element.
### 2. Its non-NA elements must be >= 1 and <= 'length(old_dim)'.
### 3. Its non-NA elements must be strictly sorted in ascending order.
### 4. 'old_dim[new2old]' must be identical to 'new_dim' after replacing
###    all NA's in it with 1's. That is, it must satisfy:
###        new_dim2 <- old_dim[new2old]
###        new_dim2[is.na(new_dim2)] <- 1L
###        identical(new_dim, new_dim2)
.map_new_to_old_dim <- function(new_dim, old_dim, x_class)
{
    cannot_map_msg <- c(
        "Cannot map the supplied dim vector to the current dimensions of ",
        "the object. On a ", x_class, " object, the dim() setter can only ",
        "be used to drop and/or add \"ineffective dimensions\" (i.e. ",
        "dimensions equal to 1) to the object."
    )
    can_map_effective_dimensions <- function(effdim_idx1, effdim_idx2) {
        if (length(effdim_idx1) != length(effdim_idx2))
            return(FALSE)
        if (length(effdim_idx1) == 0L)
            return(TRUE)
        all(new_dim[effdim_idx1] == old_dim[effdim_idx2])
    }

    ## Get index of new and old effective dimensions.
    effdim_idx1 <- which(new_dim != 1L)
    effdim_idx2 <- which(old_dim != 1L)

    if (!can_map_effective_dimensions(effdim_idx1, effdim_idx2))
        stop(wmsg(cannot_map_msg))

    ## Map as much ineffective dimensions as possible.
    ## This prevents from returning a mapping that contains only NA's, which
    ## would otherwise happen in the rare situation where 'old_dim' (and
    ## consequently 'new_dim') contains only 1's. So for example, in the
    ## following situation:
    ##     A <- DelayedArray(array(1:24, 4:2))
    ##     A1 <- A[1, 1, 1, drop=FALSE]
    ##     dim(A1) <- c(1, 1)
    ## 'A1@seed@perm' will be set to 'c(1L, 2L)', which is ok, and not to
    ## 'c(NA_integer_, NA_integer_)', which would NOT be ok.
    map_ineffective_dimensions <- function(effdim_idx1, effdim_idx2) {
        new2old <- rep.int(NA_integer_, length(new_dim))
        idx1 <- c(effdim_idx1, length(new_dim) + 1L)
        idx2 <- c(effdim_idx2, length(old_dim) + 1L)
        i1 <- i2 <- 0L
        for (k in seq_along(idx1)) {
            prev_i1 <- i1
            i1 <- idx1[[k]]
            prev_i2 <- i2
            i2 <- idx2[[k]]
            n <- min(i1 - prev_i1, i2 - prev_i2) - 1L
            new2old[prev_i1 + seq_len(n)] <- prev_i2 + seq_len(n)
        }
        new2old
    }
    new2old <- map_ineffective_dimensions(effdim_idx1, effdim_idx2)

    ## Map **all** effective dimensions.
    new2old[effdim_idx1] <- effdim_idx2
    new2old
}

.set_DelayedArray_dim <- function(x, value)
{
    x_dim <- dim(x)
    value <- unname(S4Arrays:::normalize_dim_replacement_value(value, x_dim))
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
    function(x, value) stash_DelayedSetDimnames(x, value)
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
        S4Arrays:::show_compact_array(object)
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
        S4Arrays:::combine_array_objects(objects)
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

