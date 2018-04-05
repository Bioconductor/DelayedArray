### =========================================================================
### Common operations on DelayedArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Ops" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Compare members: ==, !=, <=, >=, <, >
### Logic members: &, |
###

.normarg_Ops_vector_arg <- function(e, x_nrow,
                                    e_what="left object",
                                    x_what="first dimension of right object",
                                    x_what2=x_what,
                                    check.only=FALSE)
{
    e_len <- length(e)
    if (e_len == x_nrow || e_len == 1L)
        return(e)
    if (e_len > x_nrow)
        stop(wmsg(e_what, " is longer than ", x_what))
    if (e_len == 0L || x_nrow %% e_len != 0L)
        stop(wmsg("length of ", e_what, " is not a divisor of ", x_what2))
    if (check.only)
        return(e)
    rep(e, length.out=x_nrow)
}

### Return a DelayedArray object of the same dimensions as 'e1'.
.DelayedArray_Ops_with_right_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e1, "DelayedArray"))
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e2))
        e2 <- as.vector(e2)
    if (!is.atomic(e2))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e2 <- .normarg_Ops_vector_arg(e2, nrow(e1),
                                  e_what="right object",
                                  x_what="first dimension of left object")
    Ridx <- if (length(e2) == 1L) integer(0) else 1L
    stash_DelayedUnaryIsoOp(e1, .Generic, Rargs=list(e2), Ridx=Ridx)
}

### Return a DelayedArray object of the same dimensions as 'e2'.
.DelayedArray_Ops_with_left_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e2, "DelayedArray"))
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e1))
        e1 <- as.vector(e1)
    if (!is.atomic(e1))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e1 <- .normarg_Ops_vector_arg(e1, nrow(e2),
                                  e_what="left object",
                                  x_what="first dimension of right object")
    Lidx <- if (length(e1) == 1L) integer(0) else 1L
    stash_DelayedUnaryIsoOp(e2, .Generic, Largs=list(e1), Lidx=Lidx)
}

### Return a DelayedArray object of the same dimensions as 'e1' and 'e2'.
.DelayedArray_Ops_COMBINE_seeds <- function(.Generic, e1, e2)
{
    if (!identical(dim(e1), dim(e2)))
        stop("non-conformable arrays")
    DelayedArray(new_DelayedVariadicIsoOp(e1@seed, e2@seed, OP=.Generic))
}

.DelayedArray_Ops <- function(.Generic, e1, e2)
{
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (identical(e1_dim, e2_dim))
        return(.DelayedArray_Ops_COMBINE_seeds(.Generic, e1, e2))
    ## Effective dimensions.
    effdim_idx1 <- which(e1_dim != 1L)
    effdim_idx2 <- which(e2_dim != 1L)
    if ((length(effdim_idx1) == 1L) == (length(effdim_idx2) == 1L))
        stop("non-conformable arrays")
    if (length(effdim_idx1) == 1L) {
        .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
    } else {
        .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
    }
}

setMethod("Ops", c("DelayedArray", "vector"),
    function(e1, e2)
        .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
)

setMethod("Ops", c("vector", "DelayedArray"),
    function(e1, e2)
        .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
)

setMethod("Ops", c("DelayedArray", "DelayedArray"),
    function(e1, e2)
        .DelayedArray_Ops(.Generic, e1, e2)
)

### Support unary operators "+" and "-".
setMethod("+", c("DelayedArray", "missing"),
    function(e1, e2) stash_DelayedUnaryIsoOp(e1, .Generic, Largs=list(0L))
)
setMethod("-", c("DelayedArray", "missing"),
    function(e1, e2) stash_DelayedUnaryIsoOp(e1, .Generic, Largs=list(0L))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmax2() and pmin2()
###
### We treat them like the binary operators of the "Ops" group generics.
###

setGeneric("pmax2", function(e1, e2) standardGeneric("pmax2"))
setGeneric("pmin2", function(e1, e2) standardGeneric("pmin2"))

### Mimicking how the "Ops" members combine the "dim", "names", and "dimnames"
### attributes of the 2 operands.
.check_and_combine_dims <- function(e1, e2)
{
    dim1 <- dim(e1)
    dim2 <- dim(e2)
    if (is.null(dim1))
        return(dim2)
    if (is.null(dim2))
        return(dim1)
    if (!identical(dim1, dim2))
        stop("non-conformable arrays")
    dim1
}

.combine_names <- function(e1, e2)
{
    len1 <- length(e1)
    len2 <- length(e2)
    names1 <- names(e1)
    if (len1 > len2)
        return(names1)
    names2 <- names(e2)
    if (len2 > len1 || is.null(names1))
        return(names2)
    names1
}

setMethod("pmax2", c("ANY", "ANY"),
    function(e1, e2)
    {
        ans_dim <- .check_and_combine_dims(e1, e2)
        ans <- pmax(e1, e2)
        if (is.null(ans_dim)) {
            names(ans) <- .combine_names(e1, e2)
        } else {
            ans <- set_dim(ans, ans_dim)
            ans <- set_dimnames(ans, combine_dimnames(list(e1, e2)))
        }
        ans
    }
)

setMethod("pmin2", c("ANY", "ANY"),
    function(e1, e2)
    {
        ans_dim <- .check_and_combine_dims(e1, e2)
        ans <- pmin(e1, e2)
        if (is.null(ans_dim)) {
            names(ans) <- .combine_names(e1, e2)
        } else {
            ans <- set_dim(ans, ans_dim)
            ans <- set_dimnames(ans, combine_dimnames(list(e1, e2)))
        }
        ans
    }
)

for (.Generic in c("pmax2", "pmin2")) {
    setMethod(.Generic, c("DelayedArray", "vector"),
        function(e1, e2)
            .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("vector", "DelayedArray"),
        function(e1, e2)
            .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("DelayedArray", "DelayedArray"),
        function(e1, e2)
            .DelayedArray_Ops(.Generic, e1, e2)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sweep()
###

### Only supports a MARGIN of length 1 for now.
### Ignores 'check.margin'.
### Works if 'FUN' is a member of the Ops group or, more generally, if 'FUN'
### works on a DelayedArray object and preserves its dimensions (e.g. pmax2()
### or pmin2() above).
setMethod("sweep", "DelayedArray",
    function(x, MARGIN, STATS, FUN="-", check.margin=TRUE, ...)
    {
        FUN <- match.fun(FUN)
        if (!identical(check.margin, TRUE))
            warning(wmsg("'check.margin' is ignored when 'x' is ",
                         "a DelayedArray object or derivative"))
        x_dim <- dim(x)
        x_ndim <- length(x_dim)
        if (!isSingleNumber(MARGIN))
            stop(wmsg("the \"sweep\" method for DelayedArray objects ",
                      "only supports a MARGIN of length 1 at the moment"))
        if (!is.integer(MARGIN))
            MARGIN <- as.integer(MARGIN)
        if (MARGIN < 1 || MARGIN > x_ndim)
            stop("invalid 'MARGIN'")

        ## Check 'STATS' length.
        ## If 'FUN' is a member of the Ops group, it will check the length
        ## of 'STATS' and possibly reject it but it will display an obscure
        ## error message (see .normarg_Ops_vector_arg() in this file). By
        ## checking the length early, we can display a more appropriate
        ## error message.
        .normarg_Ops_vector_arg(STATS, x_dim[[MARGIN]],
                                e_what="'STATS'",
                                x_what="the extent of 'dim(x)[MARGIN]'",
                                x_what2="'dim(x)[MARGIN]'",
                                check.only=TRUE)

        perm <- c(MARGIN, seq_len(x_ndim)[-MARGIN])
        x2 <- aperm(x, perm)
        ans2 <- FUN(x2, STATS, ...)
        aperm(ans2, order(perm))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various unary operators + the "Math" and "Math2" groups
###
### All these operations return a DelayedArray object of the same dimensions
### as 'x'.
###

.UNARY_OPS <- c("is.na", "is.finite", "is.infinite", "is.nan", "!",
                "tolower", "toupper")

for (.Generic in .UNARY_OPS) {
    setMethod(.Generic, "DelayedArray",
        function(x) stash_DelayedUnaryIsoOp(x, .Generic)
    )
}

setMethod("lengths", "DelayedArray",
    function(x, use.names=TRUE)
        stash_DelayedUnaryIsoOp(x, "lengths",
            Rargs=list(use.names=use.names))
)

setMethod("nchar", "DelayedArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
        stash_DelayedUnaryIsoOp(x, "nchar",
            Rargs=list(type=type, allowNA=allowNA, keepNA=keepNA))
)

setMethod("Math", "DelayedArray",
    function(x) stash_DelayedUnaryIsoOp(x, .Generic)
)

.DelayedArray_Math2 <- function(.Generic, x, digits)
{
    stopifnot(is(x, "DelayedArray"))
    if (!isSingleNumberOrNA(digits))
        stop(wmsg("'digits' must be a single numeric"))
    if (!is.integer(digits))
        digits <- as.integer(digits)
    stash_DelayedUnaryIsoOp(x, .Generic, Rargs=list(digits=digits))
}

### Note that round() and signif() don't use the same default for 'digits'.
setMethod("round", "DelayedArray",
    function(x, digits=0) .DelayedArray_Math2("round", x, digits)
)
setMethod("signif", "DelayedArray",
    function(x, digits=6) .DelayedArray_Math2("signif", x, digits)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### grepl(), sub(), gsub()
###

setMethod("grepl", c(x="DelayedArray"),
    function(pattern, x,
             ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
        stash_DelayedUnaryIsoOp(x, "grepl",
            Largs=list(pattern=pattern),
            Rargs=list(ignore.case=ignore.case, perl=perl,
                       fixed=fixed, useBytes=useBytes))
)

setMethod("sub", c(x="DelayedArray"),
    function(pattern, replacement, x,
             ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
        stash_DelayedUnaryIsoOp(x, "sub",
            Largs=list(pattern=pattern, replacement=replacement),
            Rargs=list(ignore.case=ignore.case, perl=perl,
                       fixed=fixed, useBytes=useBytes))
)

setMethod("gsub", c(x="DelayedArray"),
    function(pattern, replacement, x,
             ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
        stash_DelayedUnaryIsoOp(x, "gsub",
            Largs=list(pattern=pattern, replacement=replacement),
            Rargs=list(ignore.case=ignore.case, perl=perl,
                       fixed=fixed, useBytes=useBytes))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

### Used in unit tests!
.DelayedArray_block_anyNA <- function(x, recursive=FALSE)
{
    APPLY <- anyNA
    COMBINE <- function(b, block, init, reduced) { init || reduced }
    init <- FALSE
    BREAKIF <- identity

    block_APPLY_and_COMBINE(x, APPLY, COMBINE, init, BREAKIF)
}

setMethod("anyNA", "DelayedArray", .DelayedArray_block_anyNA)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which()
###

### Used in unit tests!
.DelayedArray_block_which <- function(x, arr.ind=FALSE, useNames=TRUE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop("'arr.ind' must be TRUE or FALSE")
    if (!isTRUEorFALSE(useNames))
        stop("'useNames' must be TRUE or FALSE")
    APPLY <- base::which
    COMBINE <- function(b, block, init, reduced) {
        if (length(reduced) != 0L) {
            reduced <- reduced + init[["offset"]]
            part_number <- sprintf("%010d", b)
            init[[part_number]] <- reduced
        }
        init[["offset"]] <- init[["offset"]] + length(block)
        init
    }
    offset <- 0L
    ## If 'x' is a "long array" (i.e. longer than 2^31), we use an offset of
    ## type double to avoid integer overflow.
    x_len <- length(x)
    if (is.double(x_len))
        offset <- as.double(offset)
    init <- new.env(parent=emptyenv())
    init[["offset"]] <- offset

    init <- block_APPLY_and_COMBINE(x, APPLY, COMBINE, init)
    stopifnot(identical(x_len, init[["offset"]]))  # sanity check
    rm(list="offset", envir=init)

    if (length(init) == 0L) {
        ans <- if (is.integer(x_len)) integer(0) else numeric(0)
    } else {
        ans <- unlist(as.list(init, sorted=TRUE),
                      recursive=FALSE, use.names=FALSE)
    }
    if (arr.ind)
        ans <- arrayInd(ans, dim(x), dimnames(x), useNames=useNames)
    ans
}

setMethod("which", "DelayedArray", .DelayedArray_block_which)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max, min, range, sum, prod, any, all
###

.collect_objects <- function(x, ...)
{
    if (missing(x)) {
        objects <- unname(list(...))
    } else {
        objects <- unname(list(x, ...))
    }
    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    is_array_like <- function(x) is(x, "DelayedArray") || is.array(x)
    if (!all(vapply(objects, is_array_like, logical(1))))
        stop("the supplied objects must be array-like objects (or NULLs)")
    objects
}

### Used in unit tests!
.DelayedArray_block_Summary <- function(.Generic, x, ..., na.rm=FALSE)
{
    objects <- .collect_objects(x, ...)

    GENERIC <- match.fun(.Generic)
    APPLY <- function(block) {
        ## We get a warning if 'block' is empty (which can't happen, blocks
        ## can't be empty) or if 'na.rm' is TRUE and 'block' contains only
        ## NA's or NaN's.
        reduced <- tryCatch(GENERIC(block, na.rm=na.rm), warning=identity)
        if (is(reduced, "warning"))
            return(NULL)
        reduced
    }
    COMBINE <- function(b, block, init, reduced) {
        if (is.null(init) && is.null(reduced))
            return(NULL)
        GENERIC(init, reduced)
    }
    init <- NULL
    BREAKIF <- function(init) {
        if (is.null(init))
            return(FALSE)
        switch(.Generic,
            max=         is.na(init) || init == Inf,
            min=         is.na(init) || init == -Inf,
            range=       is.na(init[[1L]]) || all(init == c(-Inf, Inf)),
            sum=, prod=  is.na(init),
            any=         identical(init, TRUE),
            all=         identical(init, FALSE),
            FALSE)  # fallback (actually not needed)
    }

    for (x in objects)
        init <- block_APPLY_and_COMBINE(x, APPLY, COMBINE, init, BREAKIF)
    if (is.null(init))
        init <- GENERIC()
    init
}

setMethod("Summary", "DelayedArray",
    function(x, ..., na.rm=FALSE)
        .DelayedArray_block_Summary(.Generic, x, ..., na.rm=na.rm)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.DelayedArray_block_mean <- function(x, trim=0, na.rm=FALSE)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for DelayedArray objects ",
             "does not support the 'trim' argument yet")

    APPLY <- function(block) {
        tmp <- as.vector(block, mode="numeric")
        block_sum <- sum(tmp, na.rm=na.rm)
        block_nval <- length(tmp)
        if (na.rm)
            block_nval <- block_nval - sum(is.na(tmp))
        c(block_sum, block_nval)
    }
    COMBINE <- function(b, block, init, reduced) { init + reduced }
    init <- numeric(2)  # sum and nval
    BREAKIF <- function(init) is.na(init[[1L]])

    ans <- block_APPLY_and_COMBINE(x, APPLY, COMBINE, init, BREAKIF)
    ans[[1L]] / ans[[2L]]
}

### S3/S4 combo for mean.DelayedArray
mean.DelayedArray <- function(x, trim=0, na.rm=FALSE, ...)
    .DelayedArray_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "DelayedArray", .DelayedArray_block_mean)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### apply()
###

setGeneric("apply", signature="X")

.simplify_apply_answer <- function(ans)
{
    if (!all(vapply(ans, is.atomic, logical(1), USE.NAMES=FALSE)))
        return(ans)  # won't simplify

    ans_lens <- lengths(ans, use.names=FALSE)
    mat_nrow <- ans_lens[[1L]]
    if (!all(ans_lens == mat_nrow))
        return(ans)  # won't simplify

    mat_data <- unlist(unname(ans))
    if (mat_nrow == 0L)
        return(mat_data)  # zero-length atomic vector

    mat_colnames <- names(ans)
    if (mat_nrow == 1L)
        return(setNames(mat_data, mat_colnames))  # atomic vector parallel
                                                  # to 'ans'

    ## Simplify as matrix.
    mat_data_names <- names(mat_data)  # comes from the 'ans' inner names
    if (is.null(mat_data_names)) {
        mat_rownames <- NULL
    } else {
        mat_rownames <- head(mat_data_names, n=mat_nrow)
        if (!all(mat_data_names == mat_rownames))
            mat_rownames <- NULL
    }
    if (is.null(mat_rownames) && is.null(mat_colnames)) {
        mat_dimnames <- NULL
    } else {
        mat_dimnames <- list(mat_rownames, mat_colnames)
    }
    matrix(mat_data, ncol=length(ans), dimnames=mat_dimnames)
}

### MARGIN must be a single integer.
.DelayedArray_apply <- function(X, MARGIN, FUN, ...)
{
    FUN <- match.fun(FUN)
    X_dim <- dim(X)
    if (!isSingleNumber(MARGIN))
        stop("'MARGIN' must be a single integer")
    if (!is.integer(MARGIN))
        MARGIN <- as.integer(MARGIN)
    if (MARGIN < 1L || MARGIN > length(X_dim))
        stop("'MARGIN' must be >= 1 and <= length(dim(X))")

    if (X_dim[[MARGIN]] == 0L) {
        ## base::apply seems to be doing something like that!
        ans <- FUN(X, ...)
        return(as.vector(ans[0L]))
    }

    ## TODO: Try using sapply() instead of lapply(). Maybe we're lucky
    ## and it achieves the kind of simplification that we're doing with
    ## .simplify_apply_answer() so we can get rid of .simplify_apply_answer().
    ans_names <-  dimnames(X)[[MARGIN]]
    ans <- lapply(setNames(seq_len(X_dim[[MARGIN]]), ans_names),
        function(i) {
            Nindex <- vector(mode="list", length=length(X_dim))
            Nindex[[MARGIN]] <- i
            slice <- subset_by_Nindex(X, Nindex, drop=FALSE)
            slice <- set_dim(slice, dim(slice)[-MARGIN])
            FUN(slice, ...)
        })

    ## Try to simplify the answer.
    .simplify_apply_answer(ans)
}

setMethod("apply", "DelayedArray", .DelayedArray_apply)

