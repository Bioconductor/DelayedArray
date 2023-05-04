### =========================================================================
### DelayedUnaryIsoOpStack objects
### -------------------------------------------------------------------------
###
### Representation of delayed unary isometric operations stacked (a.k.a.
### piped) together.
### That is:
###
###     out <- a |> OP1 |> OP2 |> ... |> OPk
###
### where:
###   - OP1, OP2, ..., OPk are isometric array transformations i.e.
###     operations that return an array with the same dimensions as
###     the input array,
###   - 'a' is the input array,
###   - the output ('out') is an array of same dimensions as 'a'.
###
### In addition, each operation in the pipe must satisfy the property that
### each value in the output array must be determined **solely** by the
### corresponding value in the input array. In other words:
###
###     OP(a)[i_1, i_2, ..., i_n]
###
### must be equal to:
###
###     OP(a[i_1, i_2, ..., i_n])
###
### for any valid multidimensional index (i_1, i_2, ..., i_n).
###
### We refer to this property as the "locality principle".
###
### Concrete examples:
###
### 1. Things like is.na(), is.finite(), logical negation (!), nchar(),
###    tolower().
###
### 2. Most functions in the Math and Math2 groups e.g. log(), sqrt(), abs(),
###    ceiling(), round(), etc...
###    Notable exceptions are the cum*() functions (cummin(), cummax(),
###    cumsum(), and cumprod()): they don't satisfy the "locality principle".
###
### 3. Operations in the Ops group when one operand is an array and the
###    other a scalar e.g. 'a + 10', '2 ^ a', 'a <= 1', etc...
###

setClass("DelayedUnaryIsoOpStack",
    contains="DelayedUnaryIsoOp",
    representation(
        OPS="list"  # The functions to apply to the input i.e. to the
                    # incoming array-like object. For example log
                    # or function(x) log(x + 1). It should act as an
                    # isomorphism i.e. always output an array-like
                    # object **parallel** to the input (i.e. with the
                    # same dimensions as the input).
    ),
    prototype(
        OPS=list()
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

new_DelayedUnaryIsoOpStack <- function(seed=new("array"), OPS=list(),
                                       check.op=FALSE)
{
    seed_dim <- dim(seed)
    if (length(seed_dim) == 0L)
        stop(wmsg("'seed' must have dimensions"))

    if (!is.list(OPS))
        stop(wmsg("'OPS' must be a list"))
    OPS <- lapply(OPS, match.fun)

    ans <- new2("DelayedUnaryIsoOpStack", seed=seed, OPS=OPS)
    if (check.op) {
        ## We quickly test the validity of the operation by calling type()
        ## on the returned object. This will fail if the operation cannot
        ## be applied e.g. if the user does something like:
        ##   M <- DelayedArray(matrix(character(12), ncol=3))
        ##   M2 <- log(M)
        ## The test is cheap and type() will be called anyway by show()
        ## later when the user tries to display M2. Better fail early than
        ## late!
        type(ans)  # we ignore the returned value
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedUnaryIsoOpStack

.DelayedUnaryIsoOpStack_summary <- function(object)
{
    sprintf("Stack of %d unary iso op(s)", length(object@OPS))
}

summary.DelayedUnaryIsoOpStack <-
    function(object, ...) .DelayedUnaryIsoOpStack_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOpStack", summary.DelayedUnaryIsoOpStack)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###
### We inherit the default dim() and dimnames() methods defined for
### DelayedUnaryIsoOp derivatives, but overwite their extract_array() method.

setMethod("extract_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        a <- extract_array(x@seed, index)
        a_dim <- dim(a)
        for (OP in x@OPS) {
            a <- OP(a)
            ## Some operations (e.g. dnorm()) don't propagate the "dim"
            ## attribute if the input array is empty.
            a <- S4Arrays:::set_or_check_dim(a, a_dim)
        }
        a
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Propagation of sparsity
###

setMethod("is_sparse", "DelayedUnaryIsoOpStack",
    function(x)
    {
        if (!is_sparse(x@seed))
            return(FALSE)
        ## Structural sparsity will be propagated if the operations in
        ## x@OPS preserve the zeros. To find out whether zeros are preserved
        ## or not, we replace the current seed with an array of one "zero",
        ## that is, with an ordinary array of the same number of dimensions
        ## and type as the seed, but with a single "zero" element. Then we
        ## apply the operations in x@OPS to it and see whether the zero was
        ## preserved or not.
        seed_ndim <- length(dim(x@seed))
        x@seed <- S4Arrays:::make_one_zero_array(type(x@seed), seed_ndim)
        ## Same as 'as.array(x)' but doesn't try to propagate the dimnames.
        a0 <- extract_array(x, vector("list", length=seed_ndim))
        S4Arrays:::is_filled_with_zeros(a0)
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "OLD_extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("OLD_extract_sparse_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        ## Assuming that the caller respected "OLD_extract_sparse_array()
        ## Terms of Use" (see SparseArraySeed-class.R), 'is_sparse(x)'
        ## should be TRUE so we can assume that the operations in x@OPS
        ## preserve the zeros and thus only need to apply them to the
        ## nonzero data.
        sas <- OLD_extract_sparse_array(x@seed, index)
        sas_nzdata <- sas@nzdata
        for (OP in x@OPS)
            sas_nzdata <- OP(sas_nzdata)
        sas@nzdata <- sas_nzdata
        sas
    }
)

