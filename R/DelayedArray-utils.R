### =========================================================================
### Common operations on DelayedArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Binding
###
### We only support binding DelayedArray objects along the rows or the cols
### at the moment. No binding along an arbitrary dimension yet! (i.e. no
### "abind" method yet)
###

.DelayedArray_arbind <- function(...)
{
    objects <- list(...)
    stash_DelayedAbind(objects[[1L]], objects[-1L], along=1L)
}

.DelayedArray_acbind <- function(...)
{
    objects <- list(...)
    stash_DelayedAbind(objects[[1L]], objects[-1L], along=2L)
}

setMethod("arbind", "DelayedArray", .DelayedArray_arbind)
setMethod("acbind", "DelayedArray", .DelayedArray_acbind)

### Argument 'deparse.level' is ignored.
setMethod("rbind", "DelayedArray", .DelayedArray_arbind)
setMethod("cbind", "DelayedArray", .DelayedArray_acbind)

### Arguments 'use.names', 'ignore.mcols', and 'check' are ignored.
setMethod("bindROWS", "DelayedArray",
    function(x, objects=list(), use.names=TRUE, ignore.mcols=FALSE, check=TRUE)
        stash_DelayedAbind(x, objects, along=1L)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Ops" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Compare members: ==, !=, <=, >=, <, >
### Logic members: &, |
###

### Return TRUE if 'length(e)' is 1 or equal to 'nrow(x)', FALSE if it's a
### divisor of 'nrow(x)', and an error otherwise.
.check_Ops_vector_arg_length <- function(e, x_nrow,
    e_what="the left operand",
    x_what="the first dimension of the right operand",
    x_what2=x_what)
{
    e_len <- length(e)
    if (e_len == x_nrow || e_len == 1L)
        return(TRUE)
    if (is(e, "DelayedArray")) {
        e_what <- paste0("the length of ", e_what)
    } else {
        e_what <- paste0("when ", e_what, " is not a DelayedArray ",
                         "object (or derivative), its length")
    }
    if (e_len > x_nrow)
        stop(wmsg(e_what, " (", e_len, ") cannot be greater ",
                        "than ", x_what, " (", x_nrow, ")"))
    if (e_len == 0L || x_nrow %% e_len != 0L)
        stop(wmsg(e_what, " (", e_len, ") must be a divisor ",
                          "of ", x_what2, " (", x_nrow, ")"))
    FALSE
}

.normarg_Ops_vector_arg <- function(e, x_nrow,
    e_what="the left operand",
    x_what="the first dimension of the right operand")
{
    ok <- .check_Ops_vector_arg_length(e, x_nrow, e_what=e_what, x_what=x_what)
    if (!is.vector(e))
        e <- as.vector(e)  # Will realize 'e' if it's a DelayedArray object.
                           # It's important to keep the check on the length
                           # (and to fail) **before** this potentially
                           # expensive realization (e.g. if 'length(e)' is
                           # very big and greater than 'x_nrow').
    if (!ok)
        e <- rep(e, length.out=x_nrow)
    e
}

### Return a DelayedArray object of the same dimensions as 'e1'.
.DelayedArray_Ops_with_right_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e1, "DelayedArray"))
    e2 <- .normarg_Ops_vector_arg(e2, nrow(e1),
                   e_what="the right operand",
                   x_what="the first dimension of the left operand")
    if (length(e2) == 1L) {
        stash_DelayedUnaryIsoOpStack(e1,
            function(a) match.fun(.Generic)(a, e2))
    } else {
        stash_DelayedUnaryIsoOpWithArgs(e1,
            .Generic, Rargs=list(e2), Ralong=1L)
    }
}

### Return a DelayedArray object of the same dimensions as 'e2'.
.DelayedArray_Ops_with_left_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e2, "DelayedArray"))
    e1 <- .normarg_Ops_vector_arg(e1, nrow(e2),
                   e_what="the left operand",
                   x_what="the first dimension of the right operand")
    if (length(e1) == 1L) {
        stash_DelayedUnaryIsoOpStack(e2,
            function(a) match.fun(.Generic)(e1, a))
    } else {
        stash_DelayedUnaryIsoOpWithArgs(e2,
            .Generic, Largs=list(e1), Lalong=1L)
    }
}

### 'e1' and 'e2' must be DelayedArray objects. At least one of them must
### be of length 1 (i.e. have all its dimensions equal to 1) so it can be
### coerced to an ordinary vector of length 1 with as.vector(). Note that
### this coercion triggers realization.
.DelayedArray_Ops_with_a_length_one_arg <- function(.Generic, e1, e2)
{
    stopifnot(is(e1, "DelayedArray"))
    stopifnot(is(e2, "DelayedArray"))
    e1_len <- length(e1)
    e2_len <- length(e2)
    if (e1_len == 1L && e2_len == 1L) {
        ## The object with most dimensions "wins".
        e1_ndim <- length(dim(e1))
        e2_ndim <- length(dim(e2))
        if (e1_ndim > e2_ndim) {
            ## 'e1' wins.
            e2 <- as.vector(e2)  # realization
        } else if (e2_ndim > e1_ndim) {
            ## 'e2' wins.
            e1 <- as.vector(e1)  # realization
        } else {
            ## 'dim(e1)' is identical to 'dim(e2)' ==> nobody wins.
        }
    } else if (e1_len == 1L) {
        ## 'e2' wins.
        e1 <- as.vector(e1)  # realization
    } else if (e2_len == 1L) {
        ## 'e1' wins.
        e2 <- as.vector(e2)  # realization
    } else {
        ## Should never happen.
        stop(wmsg("'e1' or 'e2' must be of length 1"))
    }
    if (is.vector(e1)) {
        ## 'e2' won.
        stash_DelayedUnaryIsoOpStack(e2,
            function(a) match.fun(.Generic)(e1, a))
    } else if (is.vector(e2)) {
        ## 'e1' won.
        stash_DelayedUnaryIsoOpStack(e1,
            function(a) match.fun(.Generic)(a, e2))
    } else {
        ## Nobody won and we know that this happened because 'dim(e1)'
        ## is identical to 'dim(e2)' so 'e1' and 'e2' are conformable.
        DelayedArray(new_DelayedNaryIsoOp(.Generic, e1@seed, e2@seed))
    }
}

.DelayedArray_Ops <- function(.Generic, e1, e2)
{
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (identical(e1_dim, e2_dim))
        return(DelayedArray(new_DelayedNaryIsoOp(.Generic, e1@seed, e2@seed)))
    ## If it has only 0 or 1 effective dimensions, 'e1' or 'e2' can
    ## be treated as a vector-like argument (of length 1 if it has 0
    ## effective dimensions).
    e1_neffdim <- sum(e1_dim != 1L)
    e2_neffdim <- sum(e2_dim != 1L)
    if (e1_neffdim >= 2L && e2_neffdim >= 2L)
        stop(wmsg("non-conformable array-like objects"))
    if (e1_neffdim == 0L || e2_neffdim == 0L)
        return(.DelayedArray_Ops_with_a_length_one_arg(.Generic, e1, e2))
    if (e1_neffdim == 1L && e2_neffdim == 1L)
        stop(wmsg("non-conformable array-like objects"))
    ## The object with most effective dimensions "wins".
    if (e1_neffdim == 1L) {
        ## 'e2' wins
        ans <- .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
    } else {
        ## 'e1' wins
        ans <- .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
    }
    ans
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
    function(e1, e2)
        stash_DelayedUnaryIsoOpStack(e1, function(a) match.fun(.Generic)(a))
)
setMethod("-", c("DelayedArray", "missing"),
    function(e1, e2)
        stash_DelayedUnaryIsoOpStack(e1, function(a) match.fun(.Generic)(a))
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
            ans <- set_dimnames(ans, get_first_non_NULL_dimnames(list(e1, e2)))
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
            ans <- set_dimnames(ans, get_first_non_NULL_dimnames(list(e1, e2)))
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
        .check_Ops_vector_arg_length(STATS, x_dim[[MARGIN]],
                                     e_what="'STATS'",
                                     x_what="the extent of 'dim(x)[MARGIN]'",
                                     x_what2="'dim(x)[MARGIN]'")

        perm <- c(MARGIN, seq_len(x_ndim)[-MARGIN])
        x2 <- aperm(x, perm)
        ans2 <- FUN(x2, STATS, ...)
        aperm(ans2, order(perm))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary" and "isometric" array transformations that can be
### implemented as delayed operations (this includes the "Math" and "Math2"
### groups)
###
### All these operations return a DelayedArray object of the same dimensions
### as 'x'.
###

.UNARY_OPS <- c("is.na", "is.finite", "is.infinite", "is.nan", "!",
                "tolower", "toupper")

for (.Generic in .UNARY_OPS) {
    setMethod(.Generic, "DelayedArray",
        function(x)
            stash_DelayedUnaryIsoOpStack(x, function(a) match.fun(.Generic)(a))
    )
}

setMethod("type<-", c("DelayedArray", "ANY"),
    function(x, value)
        stash_DelayedUnaryIsoOpStack(x,
            function(a) `storage.mode<-`(a, value=value))
)

setMethod("lengths", "DelayedArray",
    function(x, use.names=TRUE)
        stash_DelayedUnaryIsoOpStack(x,
            function(a) lengths(a, use.names=use.names))
)

setMethod("nchar", "DelayedArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
        stash_DelayedUnaryIsoOpStack(x,
            function(a) nchar(a, type=type, allowNA=allowNA, keepNA=keepNA))
)

setMethod("Math", "DelayedArray",
    function(x)
        stash_DelayedUnaryIsoOpStack(x, function(a) match.fun(.Generic)(a))
)

.DelayedArray_Math2 <- function(.Generic, x, digits)
{
    stopifnot(is(x, "DelayedArray"))
    if (!isSingleNumberOrNA(digits))
        stop(wmsg("'digits' must be a single numeric"))
    if (!is.integer(digits))
        digits <- as.integer(digits)
    stash_DelayedUnaryIsoOpStack(x,
        function(a) match.fun(.Generic)(a, digits=digits))
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
        stash_DelayedUnaryIsoOpStack(x,
            function(a) grepl(pattern, a,
                              ignore.case=ignore.case, perl=perl,
                              fixed=fixed, useBytes=useBytes))
)

setMethod("sub", c(x="DelayedArray"),
    function(pattern, replacement, x,
             ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
        stash_DelayedUnaryIsoOpStack(x,
            function(a) sub(pattern, replacement, a,
                            ignore.case=ignore.case, perl=perl,
                            fixed=fixed, useBytes=useBytes))
)

setMethod("gsub", c(x="DelayedArray"),
    function(pattern, replacement, x,
             ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
        stash_DelayedUnaryIsoOpStack(x,
            function(a) gsub(pattern, replacement, a,
                             ignore.case=ignore.case, perl=perl,
                             fixed=fixed, useBytes=useBytes))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

### Used in unit tests!
.BLOCK_anyNA <- function(x, recursive=FALSE, grid=NULL)
{
    FUN <- function(block, init) {anyNA(block) || init}
    init <- FALSE
    BREAKIF <- identity
    blockReduce(FUN, x, init, BREAKIF, grid=grid)
}

.anyNA_DelayedArray <- function(x, recursive=FALSE) .BLOCK_anyNA(x, recursive)
setMethod("anyNA", "DelayedArray", .anyNA_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which()
###

.which_DelayedArray <- function(x, arr.ind=FALSE, useNames=TRUE)
{
    if (!identical(useNames, TRUE))
        warning(wmsg("'useNames' is ignored when 'x' is ",
                     "a DelayedArray object or derivative"))
    BLOCK_which(x, arr.ind=arr.ind)
}

setMethod("which", "DelayedArray", .which_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### unique() and table()
###

### We only support 1D arrays.
### Semantically equivalent to 'unique(as.array(x), ...)' which, in the 1D
### case, is also equivalent to 'unique(as.vector(x), ...)'.
### Unlike unique.array(), does not support the 'MARGIN' or 'fromLast' args.
### Return an **ordinary** 1D array.
.BLOCK_unique <- function(x, incomparables=FALSE, grid=NULL)
{
    if (length(dim(x)) != 1L)
        stop(wmsg("the \"unique\" method for DelayedArray objects ",
                  "supports 1D objects only"))

    block_results <- blockApply(x, unique, incomparables=incomparables,
                                grid=grid)

    ## Combine the block results.
    unique(unlist(block_results))
}

### S3/S4 combo for unique.DelayedArray
unique.DelayedArray <- function(x, incomparables=FALSE, ...)
    .BLOCK_unique(x, incomparables=incomparables, ...)
setMethod("unique", "DelayedArray", .BLOCK_unique)

### table()

.BLOCK_table <- function(..., grid=NULL)
{
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("the \"table\" method for DelayedArray objects ",
                  "only works on a single object at the moment"))
    x <- objects[[1L]]

    block_tables <- blockApply(x, table, grid=grid)

    ## Combine the block tables.
    levels <- unlist(lapply(block_tables, names))
    storage.mode(levels) <- type(x)
    ans_names <- as.character(unique(sort(levels)))
    block_tabs <- lapply(block_tables,
        function(block_table) {
            block_tabs <- integer(length(ans_names))
            block_tabs[match(names(block_table), ans_names)] <- block_table
            block_tabs
        })
    tab <- as.integer(rowSums(matrix(unlist(block_tabs),
                                     nrow=length(ans_names),
                                     ncol=length(block_tabs))))

    ## 'tab' is a naked integer vector. We need to decorate it (see
    ## selectMethod("table", "Rle")).
    ans_dimnames <- list(ans_names)
    names(ans_dimnames) <- S4Vectors:::.list.names(...)
    ans <- array(tab, length(tab), dimnames=ans_dimnames)
    class(ans) <- "table"
    ans
}

### The table() S4 generic is defined in BiocGenerics with dispatch on the
### ellipsis (...). Unfortunately specifying 'grid' when calling table()
### breaks dispatch. For example:
###   a <- array(sample(100L, 20000L, replace=TRUE), c(20, 4, 250))
###   A <- DelayedArray(a)
###   table(A)  # ok
###   table(A, grid=blockGrid(A, 500))
###   # Error in .BLOCK_unique(x, incomparables = incomparables, ...) : 
###   #   unused argument (nmax = nmax)
### A workaround is to call .BLOCK_table():
###   DelayedArray:::.BLOCK_table(A, grid=blockGrid(A, 500))  # ok
.table_DelayedArray <- function(...) .BLOCK_table(...)
setMethod("table", "DelayedArray", .table_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max, min, range, sum, prod, any, all
###

.collect_objects <- function(...)
{
    objects <- unname(S4Vectors:::delete_NULLs(list(...)))
    is_array_like <- function(x) is(x, "Array") || is.array(x)
    if (!all(vapply(objects, is_array_like, logical(1))))
        stop("the supplied objects must be array-like objects (or NULLs)")
    objects
}

### Used in unit tests!
### An IMPORTANT RESTRICTION is that the specified grid must be compatible
### with all the objects in '...', which means that the objects in '...'
### must be conformable!
.BLOCK_Summary <- function(.Generic, x, ..., na.rm=FALSE, grid=NULL)
{
    GENERIC <- match.fun(.Generic)
    objects <- .collect_objects(x, ...)

    FUN <- function(block, init) {
        ## We get a warning if 'block' is empty (which should happen only
        ## when 'x' itself is empty, in which case blockReduce() uses a
        ## single block that has the dimensions of 'x') or if 'na.rm' is TRUE
        ## and 'block' contains only NA's or NaN's.
        ## We use tryCatch() to catch these warnings.
        reduced_block <- tryCatch(GENERIC(block, na.rm=na.rm),
                                  warning=identity)
        if (is(reduced_block, "warning") && is.null(init))
            return(NULL)
        GENERIC(reduced_block, init)
    }
    init <- NULL
    BREAKIF <- function(init) {
        if (is.null(init))
            return(FALSE)
        switch(.Generic,
            max=         is.na(init) || (na.rm && init ==  Inf),
            min=         is.na(init) || (na.rm && init == -Inf),
            range=       is.na(init[[1L]]) ||
                             (na.rm && all(init == c(-Inf, Inf))),
            sum=, prod=  is.na(init),
            any=         identical(init, TRUE),
            all=         identical(init, FALSE),
            FALSE)  # fallback (actually not needed)
    }

    for (x in objects)
        init <- blockReduce(FUN, x, init, BREAKIF, grid)
    if (is.null(init))
        init <- GENERIC()
    init
}

.Summary_DelayedArray <- function(x, ..., na.rm=FALSE)
    .BLOCK_Summary(.Generic, x, ..., na.rm=na.rm)
setMethod("Summary", "DelayedArray", .Summary_DelayedArray)

### We override the "range" method defined above via the "Summary" method
### because we want to support the 'finite' argument like S3 method
### base::range.default() does.
### An IMPORTANT RESTRICTION is that the specified grid must be compatible
### with all the objects in '...', which means that the objects in '...'
### must be conformable!
.BLOCK_range <- function(..., na.rm=FALSE, finite=FALSE, grid=NULL)
{
    objects <- .collect_objects(...)

    FUN <- function(block, init) {
        ## See .BLOCK_Summary() above for why we use tryCatch().
        reduced_block <- tryCatch(range(block, na.rm=na.rm, finite=finite),
                                  warning=identity)
        if (is(reduced_block, "warning") && is.null(init))
            return(NULL)
        range(reduced_block, init)
    }
    init <- NULL
    BREAKIF <- function(init) {
        if (is.null(init))
            return(FALSE)
        is.na(init[[1L]]) || (na.rm && all(init == c(-Inf, Inf)))
    }

    for (object in objects)
        init <- blockReduce(FUN, object, init, BREAKIF, grid)
    if (is.null(init))
        init <- range()
    init
}

### S3/S4 combo for range.DelayedArray
range.DelayedArray <- function(..., na.rm=FALSE, finite=FALSE)
    .BLOCK_range(..., na.rm=na.rm, finite=finite)
### The signature of all the members of the S4 "Summary" group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that the S4 methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument before.
setMethod("range", "DelayedArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        .BLOCK_range(x, ..., na.rm=na.rm, finite=finite)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.BLOCK_mean <- function(x, trim=0, na.rm=FALSE, grid=NULL)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for DelayedArray objects ",
             "does not support the 'trim' argument yet")

    FUN <- function(block, init) {
        tmp <- as.vector(block, mode="numeric")
        block_sum <- sum(tmp, na.rm=na.rm)
        block_nval <- length(tmp)
        if (na.rm)
            block_nval <- block_nval - sum(is.na(tmp))
        c(block_sum, block_nval) + init
    }
    init <- numeric(2)  # sum and nval
    BREAKIF <- function(init) is.na(init[[1L]])

    ans <- blockReduce(FUN, x, init, BREAKIF, grid=grid)
    ans[[1L]] / ans[[2L]]
}

### S3/S4 combo for mean.DelayedArray
mean.DelayedArray <- function(x, trim=0, na.rm=FALSE, ...)
    .BLOCK_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "DelayedArray", .BLOCK_mean)


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
.apply_DelayedArray <- function(X, MARGIN, FUN, ...)
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
            Nindex <- vector("list", length=length(X_dim))
            Nindex[[MARGIN]] <- i
            slice <- subset_by_Nindex(X, Nindex, drop=FALSE)
            slice <- set_dim(slice, dim(slice)[-MARGIN])
            FUN(slice, ...)
        })

    ## Try to simplify the answer.
    .simplify_apply_answer(ans)
}

setMethod("apply", "DelayedArray", .apply_DelayedArray)

