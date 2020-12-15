### =========================================================================
### DelayedSetDimnames objects
### -------------------------------------------------------------------------
###
### Representation of a delayed "set dimnames" operation.
###

### Used in unit tests!
.INHERIT_FROM_SEED <- -1L

setClass("DelayedSetDimnames",
    contains="DelayedUnaryIsoOp",
    representation(
        dimnames="list"  # List with one list element per dimension in
                         # the input. Each list element must be NULL,
                         # or a character vector, or special value
                         # .INHERIT_FROM_SEED
    ),
    prototype(
        dimnames=list(.INHERIT_FROM_SEED)
    )
)

.validate_DelayedSetDimnames <- function(x)
{
    seed_dim <- dim(x@seed)
    seed_ndim <- length(seed_dim)

    ## 'dimnames' slot.
    if (length(x@dimnames) != seed_ndim)
        return(paste0("'x@dimnames' must have one list element per ",
                      "dimension in 'x@seed'"))
    ok <- mapply(function(dn, d) {
                     identical(dn, .INHERIT_FROM_SEED) ||
                     is.null(dn) ||
                     is.character(dn) && length(dn) == d
                 },
                 x@dimnames, seed_dim,
                 SIMPLIFY=FALSE, USE.NAMES=FALSE)
    if (!all(unlist(ok)))
        return(paste0("each list element in 'x@dimnames' must be NULL, ",
                      "or a character vector of length the extent of ",
                      "the corresponding dimension, or special value ",
                      .INHERIT_FROM_SEED))
    TRUE
}

setValidity2("DelayedSetDimnames", .validate_DelayedSetDimnames)

new_DelayedSetDimnames <-
    function(seed=new("array"), dimnames=.INHERIT_FROM_SEED)
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    if (identical(dimnames, .INHERIT_FROM_SEED)) {
        dimnames <- rep.int(list(.INHERIT_FROM_SEED), seed_ndim)
    } else {
        dimnames <- normarg_dimnames(dimnames, seed_dim)
        seed_dimnames <- dimnames(seed)
        dimnames <- lapply(setNames(seq_len(seed_ndim), names(dimnames)),
                           function(along) {
                               dn <- dimnames[[along]]
                               seed_dn <- seed_dimnames[[along]]
                               ## Let's play nice with seeds that return
                               ## dimnames() that are not NULL or a character
                               ## vector e.g.
                               ##   library(GDSArray)
                               ##   gds <- seqExampleFileName("gds")
                               ##   ga1 <- GDSArray(gds, "genotype/data")
                               ##   sapply(dimnames(ga1), class)
                               ##   #  variant.id   sample.id   ploidy.id
                               ##   # "character" "character"   "integer"
                               if (!(is.null(seed_dn) || is.character(seed_dn)))
                                   seed_dn <- as.character(seed_dn)
                               if (identical(dn, seed_dn))
                                   return(.INHERIT_FROM_SEED)
                               dn
                           })
    }
    new2("DelayedSetDimnames", seed=seed, dimnames=dimnames)
}

setMethod("is_noop", "DelayedSetDimnames",
    function(x)
    {
        ok <- vapply(x@dimnames, identical, logical(1), .INHERIT_FROM_SEED,
                     USE.NAMES=FALSE)
        all(ok) && identical(names(x@dimnames), names(dimnames(x@seed)))
    }
)

### S3/S4 combo for summary.DelayedSetDimnames

.DelayedSetDimnames_summary <- function(object) "Set dimnames"

summary.DelayedSetDimnames <-
    function(object, ...) .DelayedSetDimnames_summary(object, ...)

setMethod("summary", "DelayedSetDimnames", summary.DelayedSetDimnames)

### Seed contract.
### We inherit the "dim" and "extract_array" default methods for
### DelayedUnaryIsoOp derivatives, and overwite their "dimnames" method.

.get_DelayedSetDimnames_dimnames <- function(x)
{
    x_dimnames <- x@dimnames
    seed_dimnames <- dimnames(x@seed)
    ans <- lapply(setNames(seq_along(x_dimnames), names(x_dimnames)),
                  function(along) {
                      dn <- x_dimnames[[along]]
                      if (identical(dn, .INHERIT_FROM_SEED))
                          dn <- seed_dimnames[[along]]
                      dn
                  })
    simplify_NULL_dimnames(ans)
}

setMethod("dimnames", "DelayedSetDimnames", .get_DelayedSetDimnames_dimnames)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Backward compatibility with DelayedArray < 0.17.6
###
### In DelayedArray 0.17.6 the DelayedDimnames class got renamed
### DelayedSetDimnames. DelayedArray objects serialized with DelayedArray <
### 0.17.6 might contain DelayedDimnames instances nested in their "seed"
### slot so we need to keep the class around for now.
###

setClass("DelayedDimnames", contains="DelayedSetDimnames")

setMethod("updateObject", "DelayedDimnames",
    function(object, ..., verbose=FALSE)
    {
        class(object) <- "DelayedSetDimnames"
        callNextMethod()
    }
)

