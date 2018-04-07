### =========================================================================
### Visualize and simplify a tree of delayed operations
### -------------------------------------------------------------------------
###


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
        .show_tree(x@seed, indent, last.child=TRUE,
                   show.node.dim=show.node.dim)
    } else {
        if (is(x, "DelayedNaryOp"))
            x <- x@seeds
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
### simplify()
###

setGeneric("simplify", function(x) standardGeneric("simplify"))

setMethod("simplify", "ANY", identity)

setMethod("simplify", "DelayedSubset",
    function(x)
    {
        x1 <- x@seed
        if (isNoOp(x))
            return(x1)
        if (is(x1, "DelayedSubset")) {
            ## SQUASH
            index <- subset_index(x1@index, x@index)
            x <- new2("DelayedSubset", seed=x1@seed, index=index)
            return(x)
        }
        if (is(x1, "DelayedAperm")) {
            ## SWAP
            index2 <- rep.int(list(NULL), length(dim(x1@seed)))
            index2[x1@perm] <- x@index
            x2 <- new2("DelayedSubset", seed=x1@seed, index=index2)
            x1@seed <- simplify(x2)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOp")) {
            ## SWAP
            Largs <- subset_args(x1@Largs, x1@Lalong, x@index)
            Rargs <- subset_args(x1@Rargs, x1@Ralong, x@index)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- BiocGenerics:::replaceSlots(x1, seed=x,
                                                  Largs=Largs,
                                                  Rargs=Rargs)
            return(x1)
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x_dimnames <- dimnames(x)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- new_DelayedDimnames(x, x_dimnames)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedAperm",
    function(x)
    {
        x1 <- x@seed
        if (isNoOp(x))
            return(x1)
        if (is(x1, "DelayedAperm")) {
            ## SQUASH + REMOVE IF NO-OP
            x1@perm <- x1@perm[x@perm]
            if (isNoOp(x1))
                return(x1@seed)
            return(x1)
        }
        if (is(x1, "DelayedUnaryIsoOp") &&
            all(is.na(x1@Lalong)) && all(is.na(x1@Ralong))) {
                ## SWAP
                x@seed <- x1@seed
                x <- simplify(x)
                x1@seed <- x
                return(x1)
        }
        if (is(x1, "DelayedDimnames")) {
            ## SWAP
            x_dimnames <- dimnames(x)
            x@seed <- x1@seed
            x <- simplify(x)
            x1 <- new_DelayedDimnames(x, x_dimnames)
            return(x1)
        }
        x
    }
)

setMethod("simplify", "DelayedUnaryIsoOp",
    function(x)
    {
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
    function(x)
    {
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

