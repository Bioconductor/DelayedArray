### =========================================================================
### linearInd()
### -------------------------------------------------------------------------
###
### Performs the reverse transformation of arrayInd().
###


.normarg_dim2 <- function(dim)
{
    if (!is.integer(dim) || S4Vectors:::anyMissingOrOutside(dim, 0L))
        stop(wmsg("'dim' must be a vector (or matrix) ",
                  "of non-negative integers with no NAs"))
    if (!is.matrix(dim))
        dim <- matrix(dim, nrow=1L)
    dim
}

### Return an integer matrix with 1 column per dimension.
normarg_aind <- function(aind, ndim, what="'aind'")
{
    if (!is.numeric(aind))
        stop(wmsg(what, " must be a numeric vector or matrix"))
    if (is.matrix(aind)) {
        if (ncol(aind) != ndim)
            stop(wmsg(what, " must have one column (or element ",
                      "if a vector) per dimension"))
    } else {
        if (is.array(aind))
            stop(wmsg(what, " must be a numeric vector or matrix"))
        if (length(aind) != ndim)
            stop(wmsg(what, " must have one element (or column ",
                      "if a matrix) per dimension"))
        aind <- matrix(aind, ncol=ndim)
    }
    if (storage.mode(aind) != "integer")
        storage.mode(aind) <- "integer"
    aind
}

### 'aind' must be a numeric vector or matrix (a vector is treated
### like a 1-row matrix).
### Return an integer vector with one element per row in 'aind'.
linearInd <- function(aind, dim)
{
    dim <- .normarg_dim2(dim)
    ndim <- ncol(dim)
    aind <- normarg_aind(aind, ndim)
    if (nrow(dim) != 1L && nrow(dim) != nrow(aind))
        stop(wmsg("when a matrix, 'dim' must have 1 row per row in 'aind'"))
    if (ndim == 0L) {
        ans <- integer(nrow(aind))
    } else {
        ans <- aind[ , ndim]
        if (ndim >= 2L) {
            for (along in (ndim-1L):1)
                ans <- (ans - 1L) * dim[ , along] + aind[ , along]
        }
    }
    names(ans) <- rownames(aind)
    ans
}

### Return a numeric vector.
normarg_ind <- function(ind, what="'ind'")
{
    if (!is.numeric(ind))
        stop(wmsg(what, " must be a numeric vector"))
    if (is.matrix(ind))
        stop(wmsg(what, " cannot be a matrix"))
    if (suppressWarnings(min(ind, na.rm=TRUE)) < 1)
        stop(wmsg(what, " must contain positive indices"))
    ind
}

### An improved version of arrayInd() that accepts a 'dim' that is a matrix
### with 1 row per element in 'ind'.
### NOT exported.
arrayInd2 <- function(ind, dim)
{
    dim <- .normarg_dim2(dim)
    ndim <- ncol(dim)
    ind <- normarg_ind(ind)
    if (nrow(dim) != 1L && nrow(dim) != length(ind))
        stop(wmsg("when a matrix, 'dim' must have 1 row per element in 'ind'"))
    if (ndim == 0L) {
        ans <- matrix(integer(0), nrow=length(ind))
    } else {
        ind0 <- ind - 1L
        d <- dim[ , 1L]
        ans <- matrix(1L + as.integer(ind0 %% d), nrow=length(ind), ncol=ndim)
    }
    rownames(ans) <- names(ind)
    if (ndim >= 2L) {
        for (along in 2:ndim) {
            ind0 <- ind0 %/% d
            d <- dim[ , along]
            ans[ , along] <- 1L + as.integer(ind0 %% d)
        }
    }
    ans
}

