### =========================================================================
### Visualize and simplify a tree of delayed operations
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### showtree()
###
### A much more condensed version of str().
###

### Avoid use of non-ASCII characters in R source code. There must be a much
### better way to do this.
.VBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x82)))
.TEE   <- rawToChar(as.raw(c(0xe2, 0x94, 0x9c)))
.ELBOW <- rawToChar(as.raw(c(0xe2, 0x94, 0x94)))
.HBAR  <- rawToChar(as.raw(c(0xe2, 0x94, 0x80)))

### 'last.child' can be NA, TRUE, or FALSE. NA means 'x' is the root of the
### tree.
.rec_showtree <- function(x, indent="", last.child=NA, show.node.dim=TRUE)
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

        if (is(x, "DelayedOp")) {
            x_as1string <- summary(x)
        } else {
            x_as1string <- sprintf("[seed] %s object", class(x))
        }
        if (show.node.dim) {
            dim_in1string <- paste0(dim(x), collapse="x")
            x_as1string <- sprintf("%s %s: %s", dim_in1string, type(x),
                                                x_as1string)
        }
        cat(indent, prefix, x_as1string, "\n", sep="")
        if (!is(x, "DelayedOp"))
            return(invisible(NULL))
    }

    ## Display children.

    if (!is.na(last.child)) {
        ## Increase indent by 3 chars.
        indent <- paste0(indent, if (last.child) " " else .VBAR, "  ")
    }
    if (is(x, "DelayedUnaryOp")) {
        .rec_showtree(x@seed, indent, last.child=TRUE,
                      show.node.dim=show.node.dim)
    } else {
        if (is(x, "DelayedNaryOp"))
            x <- x@seeds
        nseed <- length(x)
        for (i in seq_len(nseed))
            .rec_showtree(x[[i]], indent, last.child=(i==nseed),
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
        if (isNoOp(x))
            return(x1)
        if (is(x1, "DelayedSubset")) {
            ## SQUASH + REMOVE IF NO-OP
            x1 <- subset_DelayedSubset(x1, x@index)
            if (isNoOp(x1))
                return(x1@seed)
            return(x1)
        }
        if (is(x1, "DelayedAperm")) {
            ## SWAP
            x2 <- new_DelayedSubset(x1@seed)
            x2@index[x1@perm] <- x@index
            x1@seed <- simplify(x2, incremental=TRUE)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOp")) {
            ## SWAP
            Largs <- subset_args(x1@Largs, x1@Lalong, x@index)
            Rargs <- subset_args(x1@Rargs, x1@Ralong, x@index)
            x@seed <- x1@seed
            x <- simplify(x, incremental=TRUE)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x,
                                                  Largs=Largs,
                                                  Rargs=Rargs)
            return(x1)
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

setMethod("simplify", "DelayedAperm",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
        if (isNoOp(x))
            return(x1)
        if (is(x1, "DelayedAperm")) {
            ## SQUASH + REMOVE IF NO-OP
            x1@perm <- x1@perm[x@perm]
            if (isNoOp(x1))
                return(x1@seed)
            return(simplify(x1, incremental=TRUE))
        }
        if (is(x1, "DelayedUnaryIsoOp")) {
            set_Lalong_to_NA <- !(x1@Lalong %in% x@perm)
            set_Ralong_to_NA <- !(x1@Ralong %in% x@perm)
            if (all(set_Lalong_to_NA) && all(set_Ralong_to_NA)) {
                ## SWAP
                x1@Lalong[set_Lalong_to_NA] <- NA_integer_
                x1@Ralong[set_Ralong_to_NA] <- NA_integer_
                x@seed <- x1@seed
                x <- simplify(x, incremental=TRUE)
                x1@seed <- x
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

setMethod("simplify", "DelayedUnaryIsoOp",
    function(x, incremental=FALSE)
    {
        if (!.normarg_incremental(incremental))
            x@seed <- simplify(x@seed)
        x1 <- x@seed
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
        if (isNoOp(x))
            return(x1)
        if (is(x1, "DelayedDimnames")) {
            ## SQUASH + REMOVE IF NO-OP
            x <- new_DelayedDimnames(x1@seed, dimnames(x))
            if (isNoOp(x))
                return(x@seed)
            return(x)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### netSubsetAndAperm()
###
### Only supported if nseed() == 1
###

IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS <- c(
    "is not supported on a DelayedArray object with multiple seeds at the ",
    "moment. Note that you can check the number of seeds with nseed()."
)

### Remove nodes that represent unary isomorphic ops (i.e. nodes of type
### DelayedUnaryIsoOp and DelayedDimnames) from a linear tree. Raise an error
### if the tree is not linear.
.remove_iso_ops <- function(x)
{
    if (!is(x, "DelayedOp"))
        return(x)
    if (is(x, "DelayedNaryOp")) {
        ## Tree is not linear.
        stop(wmsg("netSubsetAndAperm() ",
                  IS_NOT_SUPOORTED_IF_MULTIPLE_SEEDS))
    }
    x1 <- .remove_iso_ops(x@seed)
    if (is(x, "DelayedUnaryIsoOp") || is(x, "DelayedDimnames")) {
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
        reduced <- simplify(.remove_iso_ops(x))
        if (!is(reduced, "DelayedAperm"))
            reduced <- new_DelayedAperm(reduced)
        x1 <- reduced@seed
        if (!is(x1, "DelayedSubset"))
            reduced@seed <- new_DelayedSubset(x1)
        if (as.DelayedOp)
            return(reduced)
        ans <- reduced@seed@index
        if (!isNoOp(reduced))
            attr(ans, "dimmap") <- reduced@perm
        ans
    }
)

