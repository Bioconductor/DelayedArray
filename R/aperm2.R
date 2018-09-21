### =========================================================================
### aperm2()
### -------------------------------------------------------------------------
###
### Extend base::aperm() by allowing dropping and/or adding "ineffective
### dimensions" (i.e. dimensions equal to 1). Like base::aperm(), aperm2()
### preserves the length of 'a'.
###

### Normalize 'perm' argument of extended aperm().
### Only performs a shallow check of 'perm' e.g. won't say anything if 'perm'
### contains values > length(a_dim) or duplicates. Use validate_perm() defined
### below on the value returned by normarg_perm() for a full check.
normarg_perm <- function(perm, a_dim)
{
    if (missing(perm))
        return(rev(seq_along(a_dim)))
    if (is.null(perm))
        return(seq_along(a_dim))
    if (!is.numeric(perm))
        stop(wmsg("'perm' must be an integer vector"))
    if (!is.integer(perm))
        perm <- as.integer(perm)
    perm
}

### Validate 'perm' argument of extended aperm().
### Return TRUE if the argument is valid or a single string describing why
### it's not. This design allows the function to be used in the context of
### a validity method.
validate_perm <- function(perm, a_dim)
{
    if (!is.integer(perm))
        return("'perm' must be an integer vector")
    if (length(perm) == 0L)
        return("'perm' cannot be an empty vector")
    perm0 <- perm[!is.na(perm)]
    if (length(perm0) == 0L) {
        dropped_dims <- a_dim
    } else {
        if (S4Vectors:::anyMissingOrOutside(perm0, 1L, length(a_dim)))
            return(paste0("all non-NA values in 'perm' ",
                          "must be >= 1 and <= 'length(dim(a))'"))
        if (anyDuplicated(perm0))
            return("'perm' cannot contain non-NA duplicates")
        dropped_dims <- a_dim[-perm0]
    }
    if (!all(dropped_dims == 1L))
        return("only dimensions equal to 1 can be dropped")
    TRUE
}

### NOT exported for now (but maybe it should).
aperm2 <- function(a, perm)
{
    if (!is.array(a))
        stop(wmsg("'a' must be an array"))
    a_dim <- dim(a)
    perm <- normarg_perm(perm, a_dim)
    msg <- validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    nonNA_idx <- which(!is.na(perm))
    perm0 <- perm[nonNA_idx]
    ## We drop the dimensions not present in 'perm0'. Even though the
    ## dimensions to drop are guaranteed to be equal to 1, we should not
    ## use drop() for this because this would drop **all** the dimensions
    ## equal to 1, including some of the dimensions to keep (those can also
    ## be equal to 1).
    dim(a) <- a_dim[sort(perm0)]
    ans <- base::aperm(a, perm=rank(perm0))
    ans_dim <- rep.int(1L, length(perm))
    ans_dim[nonNA_idx] <- dim(ans)
    dim(ans) <- ans_dim
    ans
}

### Like aperm2() above the various "aperm" methods implemented in the
### package extend base::aperm() by allowing dropping and/or adding
### ineffective dimensions.
setGeneric("aperm", signature="a")

