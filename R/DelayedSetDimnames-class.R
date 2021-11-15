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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

new_DelayedSetDimnames <-
    function(seed=new("array"), dimnames=.INHERIT_FROM_SEED)
{
    seed_dim <- dim(seed)
    seed_ndim <- length(seed_dim)
    if (identical(dimnames, .INHERIT_FROM_SEED)) {
        dimnames <- rep.int(list(.INHERIT_FROM_SEED), seed_ndim)
    } else {
        dimnames <- normarg_dimnames(dimnames, seed_dim)
        ## 'dimnames(seed)' can fail e.g. if 'seed' is or contains an
        ## HDF5ArraySeed object that points to a non-existing file, but we
        ## still want new_DelayedSetDimnames() to work on such seed.
        ## Our use case for this is ExperimentHub resource EH1656. This is a
        ## SummarizedExperiment object (added to ExperimentHub on 2017-10-06
        ## by the restfulSEData folks) where the assay is a very old
        ## DelayedMatrix instance (predates DelayedArray 0.4) that binds
        ## together 14 old HDF5ArraySeed instances that point to a non-existing
        ## file ('assays.h5'). When updateObject( , check=FALSE) is called on
        ## EH1656, new_DelayedSetDimnames() gets called on a seed that contains
        ## an HDF5ArraySeed object that points to a non-existing file.
        seed_dimnames <- try(dimnames(seed), silent=TRUE)
        if (inherits(seed_dimnames, "try-error"))
            seed_dimnames <- NULL
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_noop() method
###

setMethod("is_noop", "DelayedSetDimnames",
    function(x)
    {
        ok <- vapply(x@dimnames, identical, logical(1), .INHERIT_FROM_SEED,
                     USE.NAMES=FALSE)
        ## 'dimnames(x@seed)' can fail e.g. if 'x@seed' is or contains an
        ## HDF5ArraySeed object that points to a non-existing file, but we
        ## still want is_noop() to work on such DelayedSetDimnames object.
        ## See new_DelayedSetDimnames() above for our use case.
        x_seed_dimnames <- try(dimnames(x@seed), silent=TRUE)
        if (inherits(x_seed_dimnames, "try-error"))
            x_seed_dimnames <- NULL
        all(ok) && identical(names(x@dimnames), names(x_seed_dimnames))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.DelayedSetDimnames

.DelayedSetDimnames_summary <- function(object) "Set dimnames"

summary.DelayedSetDimnames <-
    function(object, ...) .DelayedSetDimnames_summary(object, ...)

setMethod("summary", "DelayedSetDimnames", summary.DelayedSetDimnames)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Seed contract
###
### We inherit the default dim() and extract_array() methods defined for
### DelayedUnaryIsoOp derivatives, but overwite their dimnames() method.

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

