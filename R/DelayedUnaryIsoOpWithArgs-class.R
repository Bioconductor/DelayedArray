### =========================================================================
### DelayedUnaryIsoOpWithArgs objects
### -------------------------------------------------------------------------
###
### Representation of a delayed unary isometric operation with vector-like
### arguments along the dimensions of the input array.
###

setClass("DelayedUnaryIsoOpWithArgs",
    contains="DelayedUnaryIsoOp",
    representation(
        OP="function",     # The function to apply to the input array. For
                           # example `+` or `/`. It should act as an
                           # isomorphism i.e. produce an output array
                           # **parallel** to the input (i.e. with the same
                           # dimensions as the input).

        Largs="list",      # Left arguments to OP i.e. arguments to place
                           # before the input array in the function call.
        Rargs="list",      # Right arguments to OP i.e. arguments to place
                           # after the input array in the function call.

        Lalong="integer",  # One integer (or NA) per left argument indicating
                           # which dimension of the input array the argument
                           # is parallel to.
        Ralong="integer"   # One integer (or NA) per right argument indicating
                           # which dimension of the input array the argument
                           # is parallel to.
    ),
    prototype(
        OP=identity
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_Lalong_or_Ralong <- function(Lalong, Largs, seed_dim)
{
    if (identical(Lalong, NA))
        return(rep.int(NA_integer_, length(Largs)))
    if (!(is.numeric(Lalong) && length(Lalong) == length(Largs)))
        stop(wmsg("'Lalong' and 'Ralong' must be integer vectors ",
                  "parallel to 'Largs' and 'Rargs', respectively"))
    if (!is.integer(Lalong))
        Lalong <- as.integer(Lalong)
    nonNA_idx <- which(!is.na(Lalong))
    nonNA_Lalong <- Lalong[nonNA_idx]
    if (S4Vectors:::anyMissingOrOutside(nonNA_Lalong, 1L, length(seed_dim)))
        stop(wmsg("all non-NA values in 'Lalong' and 'Ralong' must ",
                  "be >= 1 and <= 'length(dim(seed))'"))
    if (any(Lalong != 1L, na.rm=TRUE))
        stop(wmsg("arguments in 'Largs' and 'Rargs' can only go along ",
                  "the first dimension of the input array at the moment"))
    ok <- elementNROWS(Largs[nonNA_idx]) == seed_dim[nonNA_Lalong]
    if (!all(ok))
        stop(wmsg("some arguments in 'Largs' and/or 'Rargs' are not ",
                  "parallel to the dimension that they go along with"))
    Lalong
}

new_DelayedUnaryIsoOpWithArgs <- function(seed=new("array"),
                                          OP=identity,
                                          Largs=list(), Rargs=list(),
                                          Lalong=NA, Ralong=NA,
                                          check.op=FALSE)
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    stopifnot(is.list(Largs), is.list(Rargs))
    Lalong <- .normarg_Lalong_or_Ralong(Lalong, Largs, seed_dim)
    Ralong <- .normarg_Lalong_or_Ralong(Ralong, Rargs, seed_dim)

    OP <- match.fun(OP)

    ans <- new2("DelayedUnaryIsoOpWithArgs", seed=seed,
                                             OP=OP,
                                             Largs=Largs, Rargs=Rargs,
                                             Lalong=Lalong, Ralong=Ralong)
    if (check.op)
        type(ans)  # we ignore the returned value
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedUnaryIsoOpWithArgs

.DelayedUnaryIsoOpWithArgs_summary <- function(object) "Unary iso op with args"

summary.DelayedUnaryIsoOpWithArgs <-
    function(object, ...) .DelayedUnaryIsoOpWithArgs_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOpWithArgs",
    summary.DelayedUnaryIsoOpWithArgs
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###
### We inherit the default dim() and dimnames() methods defined for
### DelayedUnaryIsoOp derivatives, but overwite their extract_array() method.

subset_args <- function(args, along, index)
{
    subset_arg <- function(arg, MARGIN) {
        if (is.na(MARGIN))
            return(arg)
        i <- index[[MARGIN]]
        if (is.null(i))
            return(arg)
        extractROWS(arg, i)
    }
    mapply(subset_arg, args, along, SIMPLIFY=FALSE, USE.NAMES=FALSE)
}

setMethod("extract_array", "DelayedUnaryIsoOpWithArgs",
    function(x, index)
    {
        a <- extract_array(x@seed, index)

        ## Subset the left and right arguments that go along with a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        ans <- do.call(x@OP, c(Largs, list(a), Rargs))

        ## Some operations (e.g. dnorm()) don't propagate the "dim" attribute
        ## if the input array is empty.
        .set_or_check_dim(ans, dim(a))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###
### DelayedUnaryIsoOpWithArgs objects are NOT considered to propagate
### structural sparsity at the moment. However it would be nice if
### things like 'M * runif(nrow(M))' or 'M / runif(nrow(M))' propagated
### sparsity. These are simplified versions of the following use case by
### Aaron:
###   library(TENxBrainData)
###   sce <- TENxBrainData()
###   sf <- runif(ncol(sce))
###   lcounts <- log2(t(t(counts(sce))/sf) + 1)
### 'lcounts' should be considered sparse but right now it's not!
### TODO: The is_sparse() method below should be able to propagate sparsity
### of 'A * v', 'v * A', and 'A / v',  when 'length(v)' is 'nrow(A)' and
### 'v' does not contain infinite or NA or NaN values (in the multiplication
### case) and no zero or NA or NaN values (in the division case).

setMethod("is_sparse", "DelayedUnaryIsoOpWithArgs", function(x) {
    if (is_sparse(seed(x))) {
        if (identical(x@OP, `*`)) {
            for (v in c(x@Largs, x@Rargs)) {
                if (any(!is.finite(v))) {
                    return(FALSE)
                }
            }
            return(TRUE)

        } else if (identical(x@OP, `/`)) {
            for (v in c(x@Largs, x@Rargs)) {
                if (any(!is.finite(v) | v==0)) {
                    return(FALSE)
                }
            }
            return(TRUE)
        }
    }

    FALSE
})

setMethod("extract_sparse_array", "DelayedUnaryIsoOpWithArgs",
    function(x, index) {
        ## Assuming that the caller respected "extract_sparse_array() Terms
        ## of Use" (see SparseArraySeed-class.R), 'is_sparse(x)' should be
        ## TRUE so we can assume that the operation in x@OP preserves the
        ## zeroes and thus only need to apply them to the nonzero data.
        sas <- extract_sparse_array(x@seed, index)

        ## Subset the left and right arguments that go along with a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        # Expanding to match the non-zero values.
        sas_nzindex <- sas@nzindex
        nzremap <- function(arg, MARGIN) {
            extractROWS(arg, sas_nzindex[,MARGIN])
        }
        Largs <- mapply(nzremap, arg=Largs, MARGIN=x@Lalong, SIMPLIFY=FALSE)
        Rargs <- mapply(nzremap, arg=Rargs, MARGIN=x@Ralong, SIMPLIFY=FALSE)

        sas_nzdata <- sas@nzdata
        sas_nzdata <- do.call(x@OP, c(Largs, list(sas_nzdata), Rargs))
        sas@nzdata <- sas_nzdata
        sas
})

