### =========================================================================
### DelayedArray subsetting
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BLOCK_which()
###

### Used in unit tests!
### 'x' is **trusted** to be a logical array-like object.
BLOCK_which <- function(x, arr.ind=FALSE, grid=NULL)
{
    if (!isTRUEorFALSE(arr.ind))
        stop("'arr.ind' must be TRUE or FALSE")
    ## Return a numeric matrix like one returned by base::arrayInd(), that
    ## is, a matrix where each row is an n-uplet representing an array index.
    FUN <- function(block) {
        b <- currentBlockId(block)
        m <- base::which(block)
        mapToRef(rep.int(b, length(m)), m, effectiveGrid(block), linear=TRUE)
    }
    block_results <- blockApply(x, FUN, grid=grid)
    aind <- do.call(rbind, block_results)
    aind_as_list <- lapply(ncol(aind):1, function(j) aind[ , j])
    oo <- do.call(order, aind_as_list)
    ans <- aind[oo, , drop=FALSE]
    if (!arr.ind)
        ans <- linearInd(ans, dim(x))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .BLOCK_extract_vector()
###

### Linear single bracket subsetting (e.g. x[5:2]) of an array-like object.
### Return an atomic vector.
### 'x' is **trusted** to be an array-like object.
### 'i' is **trusted** to be an integer vector representing a linear index
### of valid positions in 'x'.
.BLOCK_extract_vector <- function(x, i, grid=NULL)
{
    i_len <- length(i)
    if (i_len == 0L)
        return(as.vector(extract_empty_array(x)))
    if (i_len == 1L)
        return(extract_array_element(x, i))

    ## We don't want to use blockApply() here because it would walk on the
    ## entire grid of blocks, which is not necessary. We only need to walk
    ## on the blocks touched by linear index 'i', that is, on the blocks
    ## that contain array elements located at the positions corresponding
    ## to linear index 'i'.
    grid <- .normarg_grid(grid, x)
    nblock <- length(grid)
    x_dim <- dim(x)
    majmin <- mapToGrid(arrayInd(i, x_dim), grid, linear=TRUE)
    minor_by_block <- split(majmin$minor, majmin$major)
    res <- lapply(seq_along(minor_by_block),
        function(k) {
            b <- as.integer(names(minor_by_block)[[k]])
            m <- minor_by_block[[k]]
            if (get_verbose_block_processing())
                message("Visiting block ", b, "/", nblock, " ... ",
                        appendLF=FALSE)
            ## We don't need to load the entire block if there is only 1
            ## value to extract from it.
            if (length(m) == 1L) {
                i2 <- linearInd(mapToRef(b, m, grid, linear=TRUE), x_dim)
                block_ans <- extract_array_element(x, i2)
            } else {
                block <- read_block(x, grid[[b]])
                block_ans <- block[m]
            }
            if (get_verbose_block_processing())
                message("OK")
            block_ans
    })
    unsplit(res, majmin$major)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [
###
### Multi-dimensional form (e.g. x[3:1, , 6]) is delayed.
### Linear form (e.g. x[5:2]) is block-processed.
###

.subset_DelayedArray <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(x)  # no-op
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (nsubscript == 1L && (x_ndim != 1L || drop)) {
        ## Linear single bracket subsetting e.g. x[5:2]. NOT delayed!
        ## If 'x' is mono-dimensional and 'drop' is FALSE, we switch
        ## to "multi-dimensional single bracket subsetting" which is
        ## delayed.
        i <- Nindex[[1L]]
        if (identical(x_dim, dim(i)) && type(i) == "logical") {
            i <- BLOCK_which(i)
        } else {
            i <- normalizeSingleBracketSubscript2(i, length(x))
        }
        ## From now on 'i' is an integer vector representing a linear index
        ## of valid positions in 'x'.
        return(.BLOCK_extract_vector(x, i))
    }
    if (nsubscript != x_ndim)
        stop("incorrect number of subscripts")
    ## Multi-dimensional single bracket subsetting e.g. x[3:1, , 6]. Delayed!
    ## Return an object of the same class as 'x' (endomorphism).
    ans <- stash_DelayedSubset(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "DelayedArray", .subset_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to sparse matrix (requires the Matrix package)
###
### Based on BLOCK_which().
###

.from_DelayedMatrix_to_dgCMatrix <- function(from)
{
    idx <- BLOCK_which(from != 0L)
    array_ind <- arrayInd(idx, dim(from))
    i <- array_ind[ , 1L]
    j <- array_ind[ , 2L]
    x <- from[idx]
    Matrix::sparseMatrix(i, j, x=x, dims=dim(from), dimnames=dimnames(from))
}

setAs("DelayedMatrix", "dgCMatrix", .from_DelayedMatrix_to_dgCMatrix)
setAs("DelayedMatrix", "sparseMatrix", .from_DelayedMatrix_to_dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [<-
###

.filling_error_msg <- c(
    "filling a DelayedArray object 'x' with a value (i.e. 'x[] <- value') ",
    "is supported only when 'value' is an atomic vector and 'length(value)' ",
    "is a divisor of 'nrow(x)'"
)

.subassign_error_msg <- c(
    "subassignment to a DelayedArray object 'x' (i.e. 'x[i] <- value') is ",
    "supported only when the subscript 'i' is a logical DelayedArray object ",
    "with the same dimensions as 'x' and when 'value' is a scalar (i.e. an ",
    "atomic vector of length 1)"
)

.fill_DelayedArray_with_value <- function(x, value)
{
    if (!(is.vector(value) && is.atomic(value)))
        stop(wmsg(.filling_error_msg))
    value_len <- length(value)
    if (value_len == 1L)
        return(stash_DelayedUnaryIsoOp(x, `[<-`, Rargs=list(value=value)))
    x_len <- length(x)
    if (value_len > x_len)
        stop(wmsg("'value' is longer than 'x'"))
    x_nrow <- nrow(x)
    if (x_nrow != 0L) {
        if (value_len == 0L || x_nrow %% value_len != 0L)
            stop(wmsg(.filling_error_msg))
        value <- rep(value, length.out=x_nrow)
    }
    stash_DelayedUnaryIsoOp(x, `[<-`, Rargs=list(value=value), Ralong=1L)
}

.subassign_DelayedArray <- function(x, i, j, ..., value)
{
    if (missing(x))
        stop("'x' is missing")
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(.fill_DelayedArray_with_value(x, value))
    if (nsubscript != 1L)
        stop(wmsg(.subassign_error_msg))
    i <- Nindex[[1L]]
    if (!(is(i, "DelayedArray") &&
          identical(dim(x), dim(i)) &&
          type(i) == "logical"))
        stop(wmsg(.subassign_error_msg))
    if (!(is.vector(value) && is.atomic(value) && length(value) == 1L))
        stop(wmsg(.subassign_error_msg))
    DelayedArray(new_DelayedNaryIsoOp(x@seed, i@seed,
                                      OP=`[<-`, Rargs=list(value=value)))
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
        extract_array_element(x, i)
    }
)

