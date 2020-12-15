### =========================================================================
### DelayedUnaryIsoOpStack objects
### -------------------------------------------------------------------------
###
### Representation of delayed unary isometric operations stacked (a.k.a.
### piped) together.
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

### S3/S4 combo for summary.DelayedUnaryIsoOpStack

.DelayedUnaryIsoOpStack_summary <- function(object) "Unary iso op stack"

summary.DelayedUnaryIsoOpStack <-
    function(object, ...) .DelayedUnaryIsoOpStack_summary(object, ...)

setMethod("summary", "DelayedUnaryIsoOpStack", summary.DelayedUnaryIsoOpStack)

### Seed contract.
### We inherit the "dim" and "dimnames" default methods for DelayedUnaryIsoOp
### derivatives, and overwite their "extract_array" method.

setMethod("extract_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        a <- extract_array(x@seed, index)
        a_dim <- dim(a)
        for (OP in x@OPS) {
            a <- OP(a)
            ## Some operations (e.g. dnorm()) don't propagate the "dim"
            ## attribute if the input array is empty.
            a <- .set_or_check_dim(a, a_dim)
        }
        a
    }
)

### is_sparse() and extract_sparse_array()

### Make an ordinary array of the specified type and number of dimensions,
### and with a single "zero" element. The single element is the "zero"
### associated with the specified type e.g. the empty string ("") if type
### is "character", FALSE if it's "logical", etc... More generally, the
### "zero" element is whatever 'vector(type, length=1L)' produces.
.make_array_of_one_zero <- function(type, ndim)
{
    array(vector(type, length=1L), dim=rep.int(1L, ndim))
}

setMethod("is_sparse", "DelayedUnaryIsoOpStack",
    function(x)
    {
        if (!is_sparse(x@seed))
            return(FALSE)
        ## Structural sparsity will be propagated if the operations in
        ## x@OPS preserve the zeroes. To find out whether zeroes are preserved
        ## or not, we replace the current seed with an array of one "zero",
        ## that is, with an ordinary array of the same number of dimensions
        ## and type as the seed, but with a single "zero" element. Then we
        ## apply the operations in x@OPS to it and see whether the zero was
        ## preserved or not.
        seed_ndim <- length(dim(x@seed))
        x@seed <- .make_array_of_one_zero(type(x@seed), seed_ndim)
        a0 <- extract_array(x, rep.int(list(1L), seed_ndim))
        as.vector(a0) == vector(type(a0), length=1L)
    }
)

### 'is_sparse(x)' is assumed to be TRUE and 'index' is assumed to
### not contain duplicates. See "extract_sparse_array() Terms of Use"
### in SparseArraySeed-class.R
setMethod("extract_sparse_array", "DelayedUnaryIsoOpStack",
    function(x, index)
    {
        ## Assuming that the caller respected "extract_sparse_array() Terms
        ## of Use" (see SparseArraySeed-class.R), 'is_sparse(x)' should be
        ## TRUE so we can assume that the operations in x@OPS preserve the
        ## zeroes and thus only need to apply them to the nonzero data.
        sas <- extract_sparse_array(x@seed, index)
        sas_nzdata <- sas@nzdata
        for (OP in x@OPS)
            sas_nzdata <- OP(sas_nzdata)
        sas@nzdata <- sas_nzdata
        sas
    }
)

