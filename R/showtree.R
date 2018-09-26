### =========================================================================
### Visualize and simplify a tree of delayed operations
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
        ans <- sprintf("[seed] %s object", class(x))
    }
    if (show.node.dim) {
        dim_in1string <- paste0(dim(x), collapse="x")
        sparse <- if (is_sparse(x)) ", sparse" else ""
        ans <- sprintf("%s %s%s: %s", dim_in1string, type(x), sparse, ans)
    }
    ans
}

### Avoid use of non-ASCII characters in R source code. There must be a much
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

    if (!is.list(x)) {
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
### simplify()
###

.normarg_incremental <- function(incremental)
{
    if (!isTRUEorFALSE(incremental))
        stop("'incremental' must be TRUE or FALSE")
    incremental
}

setGeneric("simplify", signature="x",
    function(x, incremental=FALSE) standardGeneric("simplify")
)

setMethod("simplify", "ANY",
    function(x, incremental=FALSE)
    {
        .normarg_incremental(incremental)
        x
    }
)

setMethod("simplify", "DelayedSubset",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (is_noop(x))
            return(x1)
        if (is(x1, "DelayedSubset")) {
            ## MERGE + REMOVE IF NO-OP
            x1 <- subset_DelayedSubset(x1, x@index)
            if (is_noop(x1))
                return(x1@seed)
            return(x1)
        }
        if (is(x1, "DelayedAperm")) {
            ## SWAP
            index2 <- project_index_on_seed(x@index, x1)
            x2 <- new_DelayedSubset(x1@seed, index2)
            x2 <- simplify(x2, incremental=TRUE)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x2, check=FALSE)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOpStack")) {
            ## SWAP
            x2 <- x
            x2@seed <- x1@seed
            x2 <- simplify(x2, incremental=TRUE)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x2, check=FALSE)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOpWithArgs")) {
            ## SWAP
            x2 <- x
            x2@seed <- x1@seed
            x2 <- simplify(x2, incremental=TRUE)
            Largs <- subset_args(x1@Largs, x1@Lalong, x@index)
            Rargs <- subset_args(x1@Rargs, x1@Ralong, x@index)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x2,
                                                  Largs=Largs,
                                                  Rargs=Rargs,
                                                  check=FALSE)
            return(x1)
        }
        if (is(x1, "DelayedSubassign")) {
            ## SWAP
            x2 <- x
            x2@seed <- x1@seed
            x2 <- simplify(x2, incremental=TRUE)
            if (is.null(dim(x1@Rvalue))) {
                ## 'x1@Rvalue' is an ordinary vector (atomic or list) of
                ## length 1
                Lindex2 <- make_Lindex(x@index, x1@Lindex)
                Rvalue2 <- x1@Rvalue
            } else {
                ## 'x1@Rvalue' is an array-like object
                Mindex <- make_Mindex(x@index, x1@Lindex, dim(x1@seed))
                Lindex2 <- get_Lindex_from_Mindex(Mindex)
                Rindex2 <- get_Rindex_from_Mindex(Mindex)
                Rvalue2 <- new_DelayedSubset(x1@Rvalue, Rindex2)
                Rvalue2 <- simplify(Rvalue2, incremental=TRUE)
            }
            ## Even though we replace all the slots of DelayedSubassign object
            ## 'x1', this is cheaper than using new() or new_DelayedSubassign()
            ## because we can avoid the cost of validation by using
            ## 'check=FALSE'.
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x2,
                                                  Lindex=Lindex2,
                                                  Rvalue=Rvalue2,
                                                  check=FALSE)
            return(x1)
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x2 <- x
            x2@seed <- x1@seed
            x2 <- simplify(x2, incremental=TRUE)
            x1 <- new_DelayedDimnames(x2, dimnames(x))
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedAperm",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (is_noop(x))
            return(x1)
        if (is(x1, "DelayedAperm")) {
            ## MERGE + REMOVE IF NO-OP
            x1@perm <- x1@perm[x@perm]
            if (is_noop(x1))
                return(x1@seed)
            return(simplify(x1, incremental=TRUE))
        }
        if (is(x1, "DelayedUnaryIsoOpStack")) {
            ## SWAP
            x@seed <- x1@seed
            x1@seed <- simplify(x, incremental=TRUE)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOpWithArgs")) {
            perm0 <- x@perm[!is.na(x@perm)]
            set_Lalong_to_NA <- !(x1@Lalong %in% perm0)
            set_Ralong_to_NA <- !(x1@Ralong %in% perm0)
            if (all(set_Lalong_to_NA) && all(set_Ralong_to_NA)) {
                ## SWAP
                x1@Lalong[set_Lalong_to_NA] <- NA_integer_
                x1@Ralong[set_Ralong_to_NA] <- NA_integer_
                x@seed <- x1@seed
                x1@seed <- simplify(x, incremental=TRUE)
                return(x1)
            }
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x_dimnames <- dimnames(x)
            x@seed <- x1@seed
            x <- simplify(x, incremental=TRUE)
            x1 <- new_DelayedDimnames(x, x_dimnames)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedUnaryIsoOpStack",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (is(x1, "DelayedUnaryIsoOpStack")) {
            ## MERGE
            x1@OPS <- c(x1@OPS, x@OPS)
            return(x1)
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x@seed <- x1@seed
            x1@seed <- simplify(x, incremental=TRUE)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedUnaryIsoOpWithArgs",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x@seed <- x1@seed
            x1@seed <- simplify(x, incremental=TRUE)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedSubassign",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental)) {
            x@seed <- simplify(x@seed)
            x@Rvalue <- simplify(x@Rvalue)
        }
        x1 <- x@seed
        if (is_noop(x))
            return(x1)
        Rvalue <- x@Rvalue
        if (is(Rvalue, "DelayedDimnames")) {
            ## REMOVE DelayedDimnames FROM Rvalue
            x@Rvalue <- Rvalue@seed
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x@seed <- x1@seed
            x1@seed <- x
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedDimnames",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (is_noop(x))
            return(x1)
        if (is(x1, "DelayedDimnames")) {
            ## MERGE + REMOVE IF NO-OP
            x <- new_DelayedDimnames(x1@seed, dimnames(x))
            if (is_noop(x))
                return(x@seed)
            return(x)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### contentIsPristine()
###

### Return FALSE if the tree contains delayed operations that modify
### the "original array values" (i.e. the values contained in the seeds).
### The value-modifying nodes are:
###   - DelayedUnaryIsoOpStack, DelayedUnaryIsoOpWithArgs, and
###     DelayedNaryIsoOp nodes;
###   - DelayedSubassign nodes that are not no-ops.
contentIsPristine <- function(x)
{
    if (!is.list(x)) {
        if (!is(x, "DelayedOp"))
            return(TRUE)
        if (is(x, "DelayedUnaryIsoOpStack") ||
            is(x, "DelayedUnaryIsoOpWithArgs") ||
            is(x, "DelayedNaryIsoOp"))
            return(FALSE)
        if (is(x, "DelayedUnaryOp")) {
            if (is(x, "DelayedSubassign") && !is_noop(x))
                return(FALSE)
            x <- list(x@seed)
        } else {
            x <- x@seeds
        }
    }
    nchildren <- length(x)
    for (i in seq_len(nchildren)) {
        if (!contentIsPristine(x[[i]]))
            return(FALSE)
    }
    TRUE
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### netSubsetAndAperm()
###
### Only supported if nseed() == 1
###

IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS <- c(
    "is not supported on a DelayedArray object with multiple seeds at the ",
    "moment. Note that you can check the number of seeds with nseed()."
)

### Remove DelayedUnaryIsoOp nodes from a linear tree.
### Raise an error if the tree is not linear.
.remove_unary_iso_ops <- function(x)
{
    if (!is(x, "DelayedOp"))
        return(x)
    if (is(x, "DelayedNaryOp")) {
        ## Tree is not linear.
        stop(wmsg("netSubsetAndAperm() ",
                  IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS))
    }
    x1 <- .remove_unary_iso_ops(x@seed)
    if (is(x, "DelayedUnaryIsoOp")) {
        x <- x1
    } else {
        x@seed <- x1
    }
    x
}

setGeneric("netSubsetAndAperm", signature="x",
    function(x, as.DelayedOp=FALSE) standardGeneric("netSubsetAndAperm")
)

setMethod("netSubsetAndAperm", "ANY",
    function(x, as.DelayedOp=FALSE)
    {
        if (!isTRUEorFALSE(as.DelayedOp))
            stop("'as.DelayedOp' must be TRUE or FALSE")
        reduced <- simplify(.remove_unary_iso_ops(x))
        if (!is(reduced, "DelayedAperm"))
            reduced <- new_DelayedAperm(reduced)
        x1 <- reduced@seed
        if (!is(x1, "DelayedSubset"))
            reduced@seed <- new_DelayedSubset(x1)
        if (as.DelayedOp)
            return(reduced)
        ans <- reduced@seed@index
        if (!is_noop(reduced))
            attr(ans, "dimmap") <- reduced@perm
        ans
    }
)

