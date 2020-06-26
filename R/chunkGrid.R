### =========================================================================
### chunkGrid()
### -------------------------------------------------------------------------
###


### For use in *Seed classes that use a slot to store the chunkdim. See for
### example the "chunkdim" slot of the HDF5ArraySeed class defined in the
### HDF5Array package.
setClassUnion("integer_OR_NULL", c("integer", "NULL"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() generic and methods
###
### chunkdim(x) must return NULL or an integer vector compatible with dim(x).
###

### Performs the same checks as .normarg_z_chunkdim() (see AutoGrid.R).
.check_chunkdim <- function(chunkdim, x_dim)
{
    ## The relationship between 'chunkdim' and 'x_dim' must be the same
    ## as the relationship between the 'spacings' and 'refdim' slots of
    ## a RegularArrayGrid object. This guarantees that
    ## 'RegularArrayGrid(dim(x), chunkdim(x))' will always work (note
    ## that this always returns a grid with at least 1 element, even
    ## when 'x' is empty).
    if (!is.integer(chunkdim))
        return(c("didn't return an integer vector (or NULL). ",
                 "chunkdim() should always return an integer vector ",
                 "or NULL."))
    if (length(chunkdim) != length(x_dim))
        return(c("returned an integer vector of length != length(dim(x))."))
    if (S4Vectors:::anyMissingOrOutside(chunkdim, 0L))
        return("returned an integer vector with negative or NA values.")
    if (!all(chunkdim <= x_dim))
        return(c("returned chunk dimensions that are not <= their ",
                 "corresponding dimension in 'x'."))
    if (any(chunkdim == 0L & x_dim != 0L))
        return(c("returned an integer vector with illegal zeros. ",
                 "chunkdim() should always return an integer vector with ",
                 "nonzero values unless the zero values correspond to ",
                 "dimensions in 'x' that are also zero."))
    if (prod(chunkdim) > .Machine$integer.max)
        return(c("returned chunk dimensions that are too big. The ",
                 "product of the chunk dimensions should always be <= ",
                 ".Machine$integer.max"))
    TRUE
}

setGeneric("chunkdim",
    function(x)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("argument to chunkdim() must be an array-like object"))
        ans <- standardGeneric("chunkdim")
        if (is.null(ans))
            return(ans)
        msg <- .check_chunkdim(ans, x_dim)
        if (!isTRUE(msg))
            stop(wmsg("The \"chunkdim\" method for ", class(x), " objects ",
                      msg, .contact_author_msg(class(x))))
        ans
    }
)

setMethod("chunkdim", "ANY", function(x) NULL)

setMethod("chunkdim", "DelayedUnaryOp", function(x) chunkdim(x@seed))

.get_DelayedSubset_chunkdim <- function(x)
{
    seed_chunkdim <- chunkdim(x@seed)
    if (is.null(seed_chunkdim))
        return(NULL)
    ok <- lapply(seq_along(seed_chunkdim),
              function(i) {seed_chunkdim[[i]] <= 1L || is.null(x@index[[i]])})
    if (!all(unlist(ok)))
        return(NULL)
    pmin(seed_chunkdim, dim(x))
}

setMethod("chunkdim", "DelayedSubset", .get_DelayedSubset_chunkdim)

.get_DelayedAperm_chunkdim <- function(x)
{
    seed_chunkdim <- chunkdim(x@seed)
    if (is.null(seed_chunkdim))
        return(NULL)
    ans <- seed_chunkdim[x@perm]
    ans[is.na(x@perm)] <- 1L
    ans
}

setMethod("chunkdim", "DelayedAperm", .get_DelayedAperm_chunkdim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkGrid() generic and methods
###
### chunkGrid(x) must return NULL or an ArrayGrid object defining a grid on
### reference array x.
###

setGeneric("chunkGrid", function(x) standardGeneric("chunkGrid"))

setMethod("chunkGrid", "ANY",
    function(x)
    {
        x_chunkdim <- chunkdim(x)
        if (is.null(x_chunkdim))
            return(NULL)
        RegularArrayGrid(dim(x), x_chunkdim)
    }
)

