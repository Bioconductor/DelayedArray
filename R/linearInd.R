### =========================================================================
### linearInd()
### -------------------------------------------------------------------------
###
### Performs the reverse transformation of arrayInd().
###

### Return an integer matrix with 1 column per dimension.
.normarg_aind <- function(aind, ndim)
{
    if (!is.numeric(aind))
        stop(wmsg("'aind' must be a numeric vector or matrix"))
    if (is.matrix(aind)) {
        if (ncol(aind) != ndim)
            stop(wmsg("'aind' must have one column (or element ",
                      "if a vector) per dimension"))
    } else {
        if (is.array(aind))
            stop(wmsg("'aind' must be a numeric vector or matrix"))
        if (length(aind) != ndim)
            stop(wmsg("'aind' must have one element (or column ",
                      "if a matrix) per dimension"))
        aind <- matrix(aind, ncol=ndim)
    }
    if (storage.mode(aind) != "integer")
        storage.mode(aind) <- "integer"
    aind
}

### 'aind' must be a numeric vector or matrix (a vector is treated
### as a 1-row matrix).
### Return an integer vector with one element per row in 'aind'.
linearInd <- function(aind, dim)
{
    if (!is.integer(dim) || S4Vectors:::anyMissingOrOutside(dim, 0L))
        stop(wmsg("'dim' must be vector of non-negative integers"))
    ndim <- length(dim)
    aind <- .normarg_aind(aind, ndim)
    if (ndim == 0L) {
        ans <- integer(nrow(aind))
    } else {
        ans <- aind[ , ndim]
        if (ndim >= 2L) {
            for (j in (ndim-1L):1)
                ans <- (ans - 1L) * dim[[j]] + aind[ , j]
        }
    }
    names(ans) <- rownames(aind)
    ans
}

