### =========================================================================
### Visualize and access the leaves of a tree of delayed operations
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### showtree()
###
### A much more condensed version of str().
###

.node_as_one_line_summary <- function(x, show.node.dim=TRUE)
{
    if (is(x, "DelayedOp")) {
        ans <- summary(x)
    } else {
        ans <- sprintf("[seed] %s object", classNameForDisplay(x))
    }
    if (show.node.dim) {
        dim_in1string <- paste0(dim(x), collapse="x")
        sparse <- if (is_sparse(x)) ", sparse" else ""
        ans <- sprintf("%s %s%s: %s", dim_in1string, type(x), sparse, ans)
    }
    ans
}

### Avoid use of non-ASCII characters in R source code. There must be a
### better way to do this.
.VBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x82)))
.TEE   <- rawToChar(as.raw(c(0xe2, 0x94, 0x9c)))
.ELBOW <- rawToChar(as.raw(c(0xe2, 0x94, 0x94)))
.HBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x80)))

### 'last.child' can be NA, TRUE, or FALSE. NA means 'x' is the root of the
### tree.
.rec_showtree <- function(x, indent="", last.child=NA,
                             prefix="", show.node.dim=TRUE)
{
    stopifnot(isSingleString(indent))
    stopifnot(is.logical(last.child), length(last.child) == 1L)

    if (!is.list(x) || is.array(x)) {
        ## Display summary line.

        if (is.na(last.child)) {
            ## No Tprefix.
            Tprefix <- ""
        } else {
            ## 3-char Tprefix
            Tprefix <- paste0(if (last.child) .ELBOW else .TEE, .HBAR, " ")
        }
        x_as1line <- .node_as_one_line_summary(x, show.node.dim=show.node.dim)
        cat(indent, Tprefix, prefix, x_as1line, "\n", sep="")
        if (!is(x, "DelayedOp"))
            return(invisible(NULL))
    }

    ## Display children.

    if (!is.na(last.child)) {
        ## Increase indent by 3 chars.
        indent <- paste0(indent,
                         if (last.child) " " else .VBAR,
                         strrep(" ", 2 + nchar(prefix)))
    }
    if (is(x, "DelayedUnaryOp")) {
        if (is(x, "DelayedSubassign") && !is.null(dim(x@Rvalue))) {
            ## 'x@Rvalue' is an array-like object.
            .rec_showtree(x@seed, indent, last.child=FALSE,
                          show.node.dim=show.node.dim)
            .rec_showtree(x@Rvalue, indent, last.child=TRUE,
                          prefix="right value: ", show.node.dim=show.node.dim)
        } else {
            .rec_showtree(x@seed, indent, last.child=TRUE,
                          show.node.dim=show.node.dim)
        }
    } else {
        if (is(x, "DelayedNaryOp"))
            x <- x@seeds
        nchildren <- length(x)
        for (i in seq_len(nchildren))
            .rec_showtree(x[[i]], indent, last.child=(i==nchildren),
                          show.node.dim=show.node.dim)
    }
}

showtree <- function(x, show.node.dim=TRUE)
{
    if (!isTRUEorFALSE(show.node.dim))
        stop("'show.node.dim' must be TRUE or FALSE")
    .rec_showtree(x, show.node.dim=show.node.dim)
}

setMethod("show", "DelayedOp", function(object) .rec_showtree(object))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nseed()
###
### Return the number of leaves in the tree of DelayedOp nodes represented
### by 'x'. Note that nseed(x) == 1 means that the tree is linear.
###

setGeneric("nseed", function(x) standardGeneric("nseed"))

setMethod("nseed", "ANY",
    function(x)
    {
        if (is(x, "DelayedUnaryOp"))
            return(nseed(x@seed))
        if (is(x, "DelayedNaryOp"))
            x <- x@seeds
        if (is.list(x) && !is.array(x)) {
            ans <- sum(vapply(x, nseed, integer(1), USE.NAMES=FALSE))
            return(ans)
        }
        return(1L)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seedApply() and modify_seeds()
###

seedApply <- function(x, FUN, ...)
{
    if (is(x, "DelayedUnaryOp"))
        return(seedApply(x@seed, FUN, ...))
    if (is(x, "DelayedNaryOp"))
        x <- x@seeds
    if (is.list(x) && !is.array(x)) {
        ans <- lapply(x, seedApply, FUN, ...)
        return(unlist(ans, recursive=FALSE, use.names=FALSE))
    }
    list(FUN(x, ...))
}

### 'FUN' must take a seed and return a seed of the same dimensions.
### Dangerous so not intended for the end user.
### Used in HDF5Array!
modify_seeds <- function(x, FUN, ...)
{
    if (is(x, "DelayedUnaryOp")) {
        x@seed <- modify_seeds(x@seed, FUN, ...)
    } else  if (is(x, "DelayedNaryOp")) {
        x@seeds <- lapply(x@seeds, modify_seeds, FUN, ...)
    } else {
        x <- FUN(x, ...)
    }
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seed() getter/setter
###
### If nseed(x) == 1 (i.e. if 'x' is a linear tree) then the seed() getter
### and setter MUST work on 'x'.
###

IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS <- c(
    "is not supported on a DelayedArray object with multiple seeds at the ",
    "moment. Note that you can check the number of seeds with nseed()."
)

setGeneric("seed", function(x) standardGeneric("seed"))

setMethod("seed", "DelayedOp",
    function(x)
    {
        if (is(x, "DelayedNaryOp")) {
            ## Tree is not linear.
            stop(wmsg("seed() ", IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS,
                      " You can use 'seedApply(x, identity)' to extract ",
                      "all the seeds as a list."))
        }
        x1 <- x@seed
        if (!is(x1, "DelayedOp"))
            return(x1)  # found the leaf seed
        x <- x1
        callGeneric()  # recursive call
    }
)

setGeneric("seed<-", signature="x",
    function(x, ..., value) standardGeneric("seed<-")
)

.normalize_seed_replacement_value <- function(value, x_seed)
{
    x_seed_class <- class(x_seed)[[1L]]
    if (!is(value, x_seed_class))
        stop(wmsg("supplied seed must be a ", x_seed_class, " object"))
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
        seed(x1) <- value  # recursive call
        x@seed <- x1
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### The path of a DelayedOp object is the path of its leaf seed. So path()
### will work on a DelayedOp object only if it works on its leaf seed.
### For example it will work if its leaf seed is an on-disk object (e.g. an
### HDF5ArraySeed object) but not if it's an in-memory object (e.g. an
### ordinary array or RleArraySeed object).
setMethod("path", "DelayedOp",
    function(object, ...)
    {
        object_seed <- try(seed(object), silent=TRUE)
        if (is(object_seed, "try-error")) {
            ## Tree is not linear.
            stop(wmsg("path() ", IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS,
                      " You can use 'seedApply(x, path)' to extract ",
                      "all the seed paths as a list."))
        }
        path(object_seed, ...)
    }
)

### The path() setter will work on a DelayedOp object only if it works on
### its leaf seed. For example it will work if its leaf seed is an on-disk
### object (e.g. an HDF5ArraySeed object) but not if it's an in-memory object
### (e.g. an ordinary array or RleArraySeed object).
setReplaceMethod("path", "DelayedOp",
    function(object, ..., value)
    {
        object_seed <- try(seed(object), silent=TRUE)
        if (is(object_seed, "try-error")) {
            ## Tree is not linear.
            stop(wmsg("path() ", IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS))
        }
        path(object_seed, ...) <- value
        seed(object) <- object_seed
        object
    }
)

