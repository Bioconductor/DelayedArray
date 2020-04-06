### =========================================================================
### DelayedArray subsetting
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BLOCK_which()
###

.Mindex_order <- function(Mindex)
{
    cols <- lapply(ncol(Mindex):1, function(j) Mindex[ , j])
    do.call(order, cols)
}

### 'x' is **trusted** to be a logical array-like object.
### Return an L-index (if 'arr.ind=FALSE') or M-index (if 'arr.ind=TRUE').
### Used in unit tests!
BLOCK_which <- function(x, arr.ind=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(arr.ind))
        stop("'arr.ind' must be TRUE or FALSE")
    FUN <- function(block, arr.ind) {
        bid <- currentBlockId(block)
        minor <- base::which(block)
        major <- rep.int(bid, length(minor))
        grid <- effectiveGrid(block)
        Mindex <- mapToRef(major, minor, grid, linear=TRUE)
        if (arr.ind)
            return(Mindex)
        Mindex2Lindex(Mindex, refdim(grid))
    }
    block_results <- blockApply(x, FUN, arr.ind, grid=grid)
    if (arr.ind) {
        Mindex <- do.call(rbind, block_results)
        oo <- .Mindex_order(Mindex)
        ans <- Mindex[oo, , drop=FALSE]
    } else {
        ans <- sort(unlist(block_results))
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .BLOCK_subset_by_Mindex()
###
### Subset an array-like object with an M-index subscript.
### Return an ordinary vector (atomic or list).
###

### 'Mindex' is **trusted** to be a 1-row matrix representing an M-index.
### Extract a **single** array element.
.extract_array_element <- function(x, Mindex)
{
    index <- as.list(Mindex)
    a <- extract_array(x, index)
    set_dim(a, NULL)
}

### 'x' is **trusted** to be an array-like object.
### 'Mindex' is **trusted** to be an integer (or numeric) matrix
### representing an M-index.
.BLOCK_subset_by_Mindex <- function(x, Mindex, grid=NULL)
{
    ans_len <- nrow(Mindex)
    if (ans_len == 0L)
        return(vector(type(x)))
    if (ans_len == 1L)
        return(.extract_array_element(x, Mindex))

    ## We don't want to use blockApply() here because it would visit all the
    ## blocks in the grid, which is not necessary. We only need to visit the
    ## blocks touched by the L-index.
    grid <- normarg_grid(grid, x)
    nblock <- length(grid)
    majmin <- mapToGrid(Mindex, grid, linear=TRUE)
    minor_by_block <- split(majmin$minor, majmin$major)
    res <- bplapply2(seq_along(minor_by_block),
        ## TODO: Not a pure function (because it refers to 'minor_by_block',
        ## 'nblock', 'grid', and 'x') so will probably fail with
        ## parallelization backends that don't use a fork (e.g. SnowParam on
        ## Windows). Test and confirm this.
        ## FIXME: The fix is to add arguments to the function so that the
        ## objects can be passed to it.
        function(k) {
            bid <- as.integer(names(minor_by_block)[[k]])
            if (get_verbose_block_processing()) {
                message("Visiting block ", bid, "/", nblock, " ... ",
                        appendLF=FALSE)
                on.exit(message("OK"))
            }
            minor <- minor_by_block[[k]]
            ## No need to load the entire block to extract a single value
            ## from it.
            if (length(minor) == 1L) {
                Mindex1 <- mapToRef(bid, minor, grid, linear=TRUE)
                block_ans <- .extract_array_element(x, Mindex1)
            } else {
                block <- read_block(x, grid[[bid]])
                block_ans <- block[minor]
            }
            block_ans
        },
        BPPARAM=getAutoBPPARAM()
    )
    unsplit(res, majmin$major)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .BLOCK_subset_by_Lindex()
###
### Subset an array-like object with an L-index subscript.
### Return an ordinary vector (atomic or list).
###

### We only accept numeric or logical subscripts at the moment.
.normarg_Lindex <- function(Lindex, x_len)
{
    if (is.numeric(Lindex))
        return(Lindex)
    if (!is.logical(Lindex))
        stop(wmsg("invalid subscript"))
    Lindex_len <- length(Lindex)
    if (Lindex_len > x_len)
        stop(wmsg("logical subscript is longer than array to subset"))
    ## Unlike base R, we don't support NAs in logical subscripts
    ## (and if we were, we would probably just treat them as FALSE
    ## to be consistent with which()).
    if (anyNA(Lindex))
        stop(wmsg("logical subscript contains NAs"))
    if (Lindex_len == x_len || Lindex_len == 0L)
        return(which(Lindex))
    if (x_len %% Lindex_len != 0L)
        stop(wmsg("length of logical subscript (", Lindex_len, ") ",
                  "must be a divisor of array length (", x_len, ")"))
    ## Doing 'Lindex <- which(rep(Lindex, length.out=x_len))' would expand
    ## the supplied logical subscript to the length of 'x' which would be
    ## too expensive if 'Lindex' is short and 'x' is big!
    n <- x_len %/% Lindex_len
    idx <- which(Lindex)
    offsets <- rep((0L:(n-1L)) * Lindex_len, each=length(idx))
    rep.int(idx, n) + offsets
}

### 'x' must be an array-like object.
### 'Lindex' must be a valid L-index, that is, a numeric vector containing
### valid linear positions in 'x'. Logical vectors are also accepted and
### turned into a valid L-index via .normarg_Lindex().
.BLOCK_subset_by_Lindex <- function(x, Lindex, grid=NULL)
{
    x_dim <- dim(x)
    stopifnot(!is.null(x_dim))
    Lindex <- .normarg_Lindex(Lindex, prod(x_dim))
    Mindex <- Lindex2Mindex(Lindex, x_dim)
    .BLOCK_subset_by_Mindex(x, Mindex, grid=grid)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_DelayedArray_by_logical_array()
###
### Subsetting DelayedArray object 'x' by a logical array-like object with
### the same dimensions as 'x' e.g. 'x[x >= 0.9]'.
### Note that we honor 'drop' in the 1D case. This diverges from base R
### where, in the 1D case, 'x[x >= 0.9]' preserves the "dim" attribute,
### even when 'drop=TRUE' (bug?).
###

### 'x' is **trusted** to be a DelayedArray object. Note that it's only
### required to be a DelayedArray object in the 1D case and when 'drop=FALSE',
### otherwise it can be any array-like object.
### 'y' is **trusted** to be a logical array-like object of the same
### dimensions as 'x'.
.subset_DelayedArray_by_logical_array <- function(x, y, drop=TRUE)
{
    ## which() will trigger block processing (via BLOCK_which()) if 'y'
    ## is a DelayedArray object.
    if (length(dim(x)) != 1L || drop) {
        ## Return an ordinary vector (atomic or list).
        Mindex <- which(y, arr.ind=TRUE)
        .BLOCK_subset_by_Mindex(x, Mindex)
    } else {
        ## Return a 1D object of the same class and type as 'x' (endomorphism).
        Lindex <- which(y)
        stash_DelayedSubset(x, list(Lindex))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_DelayedArray_by_Mindex()
###
### Subsetting DelayedArray object 'x' by an M-index.
### Note that we honor 'drop' in the 1D case. This diverges from base R
### where, in the 1D case, 'x[matrix(5:4)]' preserves the "dim" attribute,
### even when 'drop=TRUE' (bug?).
###

### 'x' must be a DelayedArray object. Note that it's only required to be
### a DelayedArray object in the 1D case and when 'drop=FALSE', otherwise
### it can be any array-like object.
### 'Mindex' must be an integer (or numeric) matrix representing an M-index.
.subset_DelayedArray_by_Mindex <- function(x, Mindex, drop=FALSE)
{
    x_dim <- dim(x)
    stopifnot(!is.null(x_dim), is.matrix(Mindex), is.numeric(Mindex))
    if (ncol(Mindex) != length(x_dim))
        stop(wmsg("when using a numeric matrix to subset an array-like ",
                  "object, the matrix must have one column per dimension ",
                  "in the array"))
    if (length(x_dim) != 1L || drop) {
        ## Return an ordinary vector (atomic or list).
        .BLOCK_subset_by_Mindex(x, Mindex)
    } else {
        ## Return a 1D object of the same class and type as 'x' (endomorphism).
        Lindex <- Mindex2Lindex(Mindex, x_dim)
        stash_DelayedSubset(x, list(Lindex))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [
###

setMethod("drop", "DelayedArray",
    function(x)
    {
        perm <- which(dim(x) != 1L)
        if (length(perm) <= 1L)
            return(as.vector(x))  # in-memory realization
        aperm(x, perm)
    }
)

.subset_DelayedArray <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop(wmsg("'x' is missing"))
    if (!isTRUEorFALSE(drop))
        stop(wmsg("'drop' must be TRUE or FALSE"))
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(x)  # no-op
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (nsubscript == 1L) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(.subset_DelayedArray_by_logical_array(x, i, drop=drop))
        if (is.matrix(i) && is.numeric(i))
            return(.subset_DelayedArray_by_Mindex(x, i, drop=drop))
        ## Linear single bracket subsetting e.g. x[5:2].
        ## If 'x' is mono-dimensional and 'drop' is FALSE, we fallback
        ## to "multi-dimensional single bracket subsetting" which is
        ## delayed.
        if (x_ndim != 1L || drop)
            return(.BLOCK_subset_by_Lindex(x, i))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    ## Multi-dimensional single bracket subsetting e.g. x[3:1, , 6].
    ## Delayed, except when 'drop=TRUE' and the result has only one dimension,
    ## in which case the result is realized as an ordinary vector (atomic or
    ## list).
    ## In the delayed case, the result is an array-like object of the same
    ## class and type as 'x' (endomorphism), and with the same number of
    ## dimensions as 'x' (if 'drop=FALSE') or possibly less dimensions (if
    ## 'drop=TRUE').
    ans <- stash_DelayedSubset(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "DelayedArray", .subset_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to sparse matrix
###

### Calling dense2sparse() on a DelayedArray object works and triggers block
### processing twice:
###   - A 1st time when which() is called: this visits all the blocks.
###   - A 2nd time for 'x[nzindex]': this visits only the blocks with
###     nonzero values.
### See dense2sparse() implementation in R/SparseArraySeed-class.R
### .BLOCK_dense2sparse() is semantically equivalent to dense2sparse() but
### it does a single pass (all blocks are visited only once) so is more
### efficient (2x faster in average).
.BLOCK_dense2sparse <- function(x, grid=NULL)
{
    FUN <- function(block, arr.ind) {
        bid <- currentBlockId(block)
        minor <- base::which(block != 0L)
        major <- rep.int(bid, length(minor))
        grid <- effectiveGrid(block)
        nzindex <- mapToRef(major, minor, grid, linear=TRUE)
        nzdata <- block[minor]
        list(nzindex, nzdata)
    }
    block_results <- blockApply(x, FUN, grid=grid)
    nzindex_list <- lapply(block_results, `[[`, 1L)
    nzdata_list <- lapply(block_results, `[[`, 2L)
    nzindex <- do.call(rbind, nzindex_list)
    nzdata <- unlist(nzdata_list, recursive=FALSE)
    SparseArraySeed(dim(x), nzindex, nzdata, check=FALSE)
}

setAs("DelayedArray", "SparseArraySeed",
    function(from) .BLOCK_dense2sparse(from)
)

.from_DelayedMatrix_to_dgCMatrix <- function(from)
{
    ans <- as(.BLOCK_dense2sparse(from), "dgCMatrix")
    ## The above does NOT propagate the dimnames at the moment (the
    ## intermediate SparseArraySeed container cannot currently store them)
    ## so we propagate them explicitly.
    from_dimnames <- dimnames(from)
    if (!is.null(from_dimnames))
        dimnames(ans) <- from_dimnames
    ans
}
setAs("DelayedMatrix", "dgCMatrix", .from_DelayedMatrix_to_dgCMatrix)
setAs("DelayedMatrix", "sparseMatrix", .from_DelayedMatrix_to_dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [<- (a.k.a. subassignment)
###
### Supported forms:
###
### - Multi-dimensional form (e.g. x[3:1, , 6] <- 0) is supported.
###
### - Linear form (i.e. x[i] <- value) is supported only when the
###   subscript 'i' is a logical DelayedArray object of the same
###   dimensions as 'x' and 'value' an ordinary vector of length 1.
###   So for example x[x <= 0.05] <- 0 is supported but x[5:2] <- 0
###   is not at the moment.
###
### - Filling (i.e. x[] <- value) is supported only when 'value' is an
###   ordinary vector with a length that is a divisor of 'nrow(x)'.
###
### All these forms of subassignment are delayed.
###


.linear_subassign_error_msg <- c(
    "linear subassignment to a DelayedArray object 'x' (i.e. 'x[i] <- value') ",
    "is only supported when the subscript 'i' is a logical DelayedArray ",
    "object of the same dimensions as 'x' and 'value' an ordinary vector ",
    "of length 1)"
)

.filling_error_msg <- c(
    "filling a DelayedArray object 'x' with a vector 'v' (i.e. 'x[] <- v') ",
    "is only supported if 'v' is an ordinary vector of length a divisor ",
    "of 'nrow(x)'"
)

.fill_DelayedArray_with_vector <- function(x, v)
{
    stopifnot(is.vector(v))
    x_len <- length(x)
    v_len <- length(v)
    if (v_len > x_len)
        stop(wmsg("right value is longer than left value"))
    x_nrow <- nrow(x)
    if (x_nrow != 0L) {
        if (v_len == 0L || x_nrow %% v_len != 0L)
            stop(wmsg(.filling_error_msg))
        v <- rep(v, length.out=x_nrow)
    }
    stash_DelayedUnaryIsoOpWithArgs(x, `[<-`,
                                    Rargs=list(value=v), Ralong=1L)
}

.subassign_DelayedArray <- function(x, i, j, ..., value)
{
    if (missing(x))
        stop("'x' is missing")
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (nsubscript == 1L && x_ndim != 1L) {
        ## Linear form.
        i <- Nindex[[1L]]
        if (!(is(i, "DelayedArray") &&
              identical(x_dim, dim(i)) &&
              type(i) == "logical"))
            stop(wmsg(.linear_subassign_error_msg))
        if (!(is.vector(value) && length(value) == 1L))
            stop(wmsg(.linear_subassign_error_msg))
        ans <- DelayedArray(new_DelayedNaryIsoOp(`[<-`, x@seed, i@seed,
                                                 Rargs=list(value=value)))
        return(ans)
    }
    ## Multi-dimensional form.
    if (nsubscript == 0L) {
        Nindex <- vector("list", length=x_ndim)
        fill_mode <- TRUE
    } else {
        if (nsubscript != x_ndim)
            stop("incorrect number of subscripts")
        fill_mode <- all(S4Vectors:::sapply_isNULL(Nindex))
    }
    if (fill_mode && is.vector(value) && length(value) != 1L)
        return(.fill_DelayedArray_with_vector(x, value))
    stash_DelayedSubassign(x, Nindex, value)
}

setReplaceMethod("[", "DelayedArray", .subassign_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [[
###
### TODO: Implement a "getArrayElement" method for Array objects (in
### extract_array.R) that is based on extract_array() e.g it just needs
### to do something like:
###
###     as.vector(extract_array(x, as.list(subscripts)))
###
### This will make multi-dimensional and linear [[ work on DelayedArray
### objects. Then remove the method below and update DelayedArray-class.Rd
###

### Only support linear subsetting at the moment.
### TODO: Support multi-dimensional subsetting e.g. x[[5, 15, 2]] or
### x[["E", 15, "b"]].
setMethod("[[", "DelayedArray",
    function(x, i, j, ...)
    {
        dots <- list(...)
        if (length(dots) > 0L)
            dots <- dots[names(dots) != "exact"]
        if (!missing(j) || length(dots) > 0L)
            stop("incorrect number of subscripts")
        if (!is.numeric(i))
            stop("invalid [[ subscript type: ", class(i)[[1L]])
        if (length(i) < 1L)
            stop("attempt to extract less than one element")
        if (length(i) > 1L)
            stop("attempt to extract more than one element")
        if (is.na(i))
            stop("NA is not a valid [[ subscript")
        Mindex <- Lindex2Mindex(i, dim(x))
        .extract_array_element(x, Mindex)[[1L]]
    }
)

