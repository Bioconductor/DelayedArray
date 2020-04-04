### =========================================================================
### Manipulation of array selections
### -------------------------------------------------------------------------
###


### Like base::arrayInd() but faster and accepts a matrix for 'dim' (with 1
### row per element in 'Lindex').
Lindex2Mindex <- function(Lindex, dim, use.names=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    .Call2("C_Lindex2Mindex", Lindex, dim, use.names, PACKAGE="DelayedArray")
}

Mindex2Lindex <- function(Mindex, dim, use.names=FALSE, as.integer=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.integer))
        stop("'as.integer' must be TRUE or FALSE")
    if (storage.mode(Mindex) != "integer")
        storage.mode(Mindex) <- "integer"
    .Call2("C_Mindex2Lindex", Mindex, dim, use.names,
                              as.integer, PACKAGE="DelayedArray")
}

### 'aind' must be an integer matrix or vector (a vector is treated like
### a 1-row matrix).
### Return a numeric vector with one element per row in 'aind'.
linearInd <- function(aind, dim)
{
    .Deprecated("Mindex2Lindex")
    Mindex2Lindex(aind, dim, use.names=TRUE)
}

