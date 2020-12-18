### =========================================================================
### DelayedUnaryIsoOpWithArgs objects
### -------------------------------------------------------------------------
###
### Representation of a delayed unary isometric operation with vector-like
### arguments going along the dimensions of the input array.
### That is:
###
###     out <- OP(L1, L2, ..., a, R1, R2, ...)
###
### where:
###   - OP is an isometric array transformation i.e. an operation that
###     returns an array with the same dimensions as the input array.
###   - 'a' is the input array,
###   - 'L1', 'L2', ..., are the left arguments,
###   - 'R1', 'R2', ..., are the right arguments,
###   - the output ('out') is an array of same dimensions as 'a'.
###
### Some of the arguments (left or right) can be vector-like arguments that
### go along the dimensions of the input array. For example if 'a' is a
### 12 x 150 x 5 array, argument 'L2' is considered to go along the 3rd
### dimension if its length is 5 and if the result of:
###
###     OP(L1, L2[k], ..., a[ , , k, drop=FALSE], R1, R2, ...)
###
### is the same as 'out[ , , k, drop=FALSE]' for any index 'k'.
###
### More generally speaking, if, say, arguments 'L2', 'L3', 'R1', and 'R2'
### go along the 3rd, 1st, 2nd, and 1st dimensions, respectively, then each
### value in the output array ('out[i, j, k]') must be determined solely by
### the corresponding values in the input array ('a[i, j, k]') and arguments
### ('L2[k]', 'L3[i]', 'R1[j]', 'R2[i]'). In other words, 'out[i, j, k]'
### must be equal to:
###
###    OP(L1, L2[k], L3[i], ..., a[i, j, k], R1[j], R2[i], ...)
###
### for any 1 <= 'i' <= 12, 1 <= 'j' <= 150, and 1 <= 'k' <= 5.
###
### We refer to this property as the "locality principle".
###
### Concrete examples:
### 1. Addition (or any operation in the Ops group) of an array 'a' and an
###    atomic vector 'v' of length 'dim(a)[[1]]':
###    - `+`(a, v):  OP is `+`, right argument goes along the 1st dimension.
###    - `<=`(a, v): OP is `<=`, right argument goes along the 1st dimension.
###    - `&`(v, a):  OP is `&`, left argument goes along the 1st dimension.
### 2. scale(x, center=v1, scale=v2): OP is `scale`, right arguments 'center'
###    and 'scale' go along the 2nd dimension.
###
### Note that if OP has no argument that goes along a dimension of
### the input array, then the delayed operation is better represented with
### a DelayedUnaryIsoOpStack object.
###

setClass("DelayedUnaryIsoOpWithArgs",
    contains="DelayedUnaryIsoOp",
    representation(
        ## 'OP' is the function to apply to the input array. For example `+`
        ## or `<=`. Must be an isometric array transformation that satisfies
        ## the "locality principle" (see above).
        OP="function",

        ## 'Largs' and 'Rargs' are the left and right arguments to 'OP()',
        ## respectively, i.e. the arguments to place before and after the
        ## input array in the function call.
        Largs="list",
        Rargs="list",

        ## 'Lalong' and 'Ralong' are integer vectors parallel to 'Largs' and
        ## 'Rargs' respectively. 'Lalong[i]' indicates which dimension of
        ## the input array the i-th left-argument ('Largs[[i]]') goes along.
        ## An NA means that the argument doesn't go along any dimension.
        Lalong="integer",
        Ralong="integer"
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

        ## Subset the left and right arguments that go along a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        ans <- do.call(x@OP, c(Largs, list(a), Rargs))

        ## Some operations (e.g. dnorm()) don't propagate the "dim" attribute
        ## if the input array is empty.
        set_or_check_dim(ans, dim(a))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

### DelayedUnaryIsoOpWithArgs object 'x' is considered to propagate sparsity
### iff the zeros in 'x@seed' are realized as zeros in 'as.array(x)'.
### For example if 'x@seed' is a 12 x 150 x 5 array and 'x@Larg[[2]]',
### 'x@Larg[[3]]', 'x@Rarg[[1]]', and 'x@Rarg[[2]]', are vectors that go
### along the 3rd, 1st, 2nd, and 1st dimensions, respectively, then by virtue
### of the "locality principle" (see at top of this file), 'x' is considered
### to propagate sparsity iff:
###
###    x@OP(x@Larg[[1]], x@Larg[[2]][k], x@Larg[[3]][i], ...,
###         x@seed[i, j, k],
###         x@Rarg[[1]][j], x@Rarg[[2]][i], ...)
###
### is a zero for any valid 3D index (i, j, k) for which 'a[i, j, k]' is a
### zero. However, performing a test like this is equivalent to computing
### the full output array which is not an option in general as it would
### defeat the purpose of using delayed operations.
###
### So we use the following simplified test instead:
###
###   1. If 'x' has arguments that go along more than one dimension, we give
###      up and declare that sparsity is not propagated.
###
###   2. Assuming that all the arguments in 'x' that go along a dimension go
###      along **the same dimension**, say, the p-th dimension, then we can
###      create a zero-filled ordinary array 'seed0' with the same number of
###      dimensions as 'x@seed' but where all the dimensions are set to 1
###      except the p-th dimension which we set to 'dim(x@seed)[[along]]'.
###      Note that 'seed0' is parallel to all the vector-like arguments that
###      go along the p-th dimension. Then if:
###
###        x@OP(x@Larg[[1]], x@Larg[[2]], ...,
###             seed0,
###             x@Rarg[[1]], x@Rarg[[2]], ...)
###
###      is an array (of the same geometry as 'seed0') filled with zeros,
###      then we know that 'x' propagates zeros.
###
### Note that this test is simple and fast BUT it can produce false negatives,
### that is, it cannot detect all the situations where sparsity is propagated.
setMethod("is_sparse", "DelayedUnaryIsoOpWithArgs",
    function(x)
    {
        if (!is_sparse(x@seed))
            return(FALSE)
        p <- setdiff(c(x@Lalong, x@Ralong), NA_integer_)
        if (length(p) >= 2L)
            return(FALSE)
        seed_ndim <- length(dim(x@seed))
        dim0 <- rep.int(1L, seed_ndim)
        if (length(p) == 1L)
            dim0[[p]] <- dim(x@seed)[[p]]
        x@seed <- make_zero_filled_array(type(x@seed), dim0)
        ## Same as 'as.array(x)' but doesn't try to propagate the dimnames.
        a0 <- extract_array(x, vector("list", length=seed_ndim))
        is_filled_with_zeros(a0)
    }
)

setMethod("extract_sparse_array", "DelayedUnaryIsoOpWithArgs",
    function(x, index)
    {
        ## Assuming that the caller respected "extract_sparse_array() Terms
        ## of Use" (see SparseArraySeed-class.R), 'is_sparse(x)' should be
        ## TRUE so we can assume that the operation in x@OP preserves the
        ## zeros and thus only need to apply them to the nonzero data.
        sas <- extract_sparse_array(x@seed, index)

        ## Subset the left and right arguments that go along a dimension.
        Largs <- subset_args(x@Largs, x@Lalong, index)
        Rargs <- subset_args(x@Rargs, x@Ralong, index)

        ## Expanding to match the non-zero values.
        sas_nzindex <- sas@nzindex
        nzremap <- function(arg, MARGIN) {
            extractROWS(arg, sas_nzindex[,MARGIN])
        }
        Largs <- mapply(nzremap, arg=Largs, MARGIN=x@Lalong, SIMPLIFY=FALSE)
        Rargs <- mapply(nzremap, arg=Rargs, MARGIN=x@Ralong, SIMPLIFY=FALSE)

        sas@nzdata <- do.call(x@OP, c(Largs, list(sas@nzdata), Rargs))
        sas
    }
)

