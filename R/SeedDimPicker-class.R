### =========================================================================
### SeedDimPicker objects
### -------------------------------------------------------------------------
###
### This class is for internal use only and is not exported.
###

setClass("SeedDimPicker",
    representation(
        seed="ANY",                # An array-like object expected to satisfy
                                   # the "seed contract" i.e. to support dim(),
                                   # dimnames(), and extract_array().

        dim_combination="integer"  # Index into dim(seed) specifying the seed
                                   # dimensions to keep.
    ),
    prototype(
        seed=new("array"),
        dim_combination=1L
    )
)

.validate_SeedDimPicker <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)
    ## 'seed' slot.
    if (seed_ndim == 0L)
        return(wmsg2("'x@seed' must have dimensions"))
    ## 'dim_combination' slot.
    if (length(x@dim_combination) == 0L)
        return(wmsg2("'x@dim_combination' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@dim_combination, 1L, seed_ndim))
        return(wmsg2("all values in 'x@dim_combination' must be >= 1 ",
                     "and <= 'seed_ndim'"))
    if (anyDuplicated(x@dim_combination))
        return(wmsg2("'x@dim_combination' cannot have duplicates"))
    if (!all(seed_dim[-x@dim_combination] == 1L))
        return(wmsg2("dimensions to drop from 'x' must be equal to 1"))
    TRUE
}

setValidity2("SeedDimPicker", .validate_SeedDimPicker)

new_SeedDimPicker <- function(seed, dim_combination=NULL)
{
    seed <- remove_pristine_DelayedArray_wrapping(seed)
    if (is.null(dim_combination))
        dim_combination <- seq_along(dim(seed))
    new2("SeedDimPicker", seed=seed, dim_combination=dim_combination)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Implement the "seed contract" i.e. dim(), dimnames(), and extract_array()
###

.get_SeedDimPicker_dim <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_dim[x@dim_combination]
}

setMethod("dim", "SeedDimPicker", .get_SeedDimPicker_dim)

.get_SeedDimPicker_dimnames <- function(x)
{
    seed_dimnames <- dimnames(x@seed)
    seed_dimnames[x@dim_combination]  # return NULL if 'seed_dimnames' is NULL
}

setMethod("dimnames", "SeedDimPicker", .get_SeedDimPicker_dimnames)

.extract_array_from_SeedDimPicker <- function(x, index)
{
    seed_dim <- dim(x@seed)
    seed_index <- rep.int(list(1L), length(seed_dim))
    seed_index[x@dim_combination] <- index
    subseed <- extract_array(x@seed, seed_index)
    dim(subseed) <- dim(subseed)[sort(x@dim_combination)]
    aperm(subseed, perm=rank(x@dim_combination))
}

setMethod("extract_array", "SeedDimPicker",
    .extract_array_from_SeedDimPicker
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###

setGeneric("aperm", signature="a")

### Unlike base::aperm() the method below supports dropping dimensions.
### If 'simplify' is TRUE, 'aperm(a)' drops the SeedDimPicker wrapping
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

.aperm.SeedDimPicker <- function(a, perm, simplify=TRUE)
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

### S3/S4 combo for aperm.SeedDimPicker
aperm.SeedDimPicker <- function(a, perm, ...)
    .aperm.SeedDimPicker(a, perm, ...)
setMethod("aperm", "SeedDimPicker", aperm.SeedDimPicker)

