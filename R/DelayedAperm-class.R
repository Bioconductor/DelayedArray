### =========================================================================
### DelayedAperm objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###

setGeneric("aperm", signature="a")

### Unlike base::aperm() the method below supports dropping dimensions.
### If 'simplify' is TRUE, 'aperm(a)' drops the DelayedAperm wrapping
### around the returned object if this wrapping represents a dim combination
### that is the identity (i.e. if the wrapped seed is semantically equivalent
### to the seed).

normarg_perm <- function(perm, a_dim)
{
    if (!is.numeric(perm))
        stop(wmsg("'perm' must be an integer vector"))
    if (!is.integer(perm))
        perm <- as.integer(perm)
    if (length(perm) == 0L)
        stop(wmsg("'perm' cannot be an empty vector"))
    if (S4Vectors:::anyMissingOrOutside(perm, 1L, length(a_dim)))
        stop(wmsg("values out of range in 'perm'"))
    if (anyDuplicated(perm))
        stop(wmsg("'perm' cannot have duplicates"))
    if (!all(a_dim[-perm] == 1L))
        stop(wmsg("dimensions to drop from 'a' must be equal to 1"))
    perm
}

.aperm.DelayedAperm <- function(a, perm, simplify=TRUE)
{
    if (!isTRUEorFALSE(simplify))
        stop(wmsg("'simplify' must be TRUE or FALSE"))
    if (missing(perm)) {
        a@dim_combination <- rev(a@dim_combination)
    } else {
        perm <- normarg_perm(perm, dim(a))
        ## At this point we know that 'perm' is a valid combination so
        ## we don't need to validate 'a' after the change below.
        a@dim_combination <- a@dim_combination[perm]
    }
    if (simplify && identical(a@dim_combination, seq_along(dim(a@seed))))
        return(a@seed)
    a
}

### S3/S4 combo for aperm.DelayedAperm
aperm.DelayedAperm <- function(a, perm, ...)
    .aperm.DelayedAperm(a, perm, ...)
setMethod("aperm", "DelayedAperm", aperm.DelayedAperm)

