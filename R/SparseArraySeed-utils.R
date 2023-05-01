### =========================================================================
### Operate natively on SparseArraySeed objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Binding
###

### Similar to S4Arrays:::simple_abind() (see R/abind.R in the S4Arrays
### package) but works on SparseArraySeed objects.
abind_SparseArraySeed_objects <- function(objects, along)
{
    stopifnot(is.list(objects))
    if (length(objects) == 0L)
        return(NULL)

    ## Check dim compatibility.
    dims <- S4Arrays:::get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Compute 'ans_dim' and 'ans_dimnames'.
    ans_dim <- S4Arrays:::combine_dims_along(dims, along)
    ans_dimnames <- S4Arrays:::combine_dimnames_along(objects, dims, along)

    ## Combine the "nzindex" slots.
    offsets <- cumsum(dims[along, -ncol(dims)])
    nzindex_list <- lapply(seq_along(objects),
        function(i) {
            object <- objects[[i]]
            nzindex <- object@nzindex
            if (i >= 2L)
                nzindex[ , along] = nzindex[ , along, drop=FALSE] +
                                    offsets[[i - 1L]]
            nzindex
	}
    )
    ans_nzindex <- do.call(rbind, nzindex_list)

    ## Combine the "nzdata" slots.
    ans_nzdata <- unlist(lapply(objects, slot, "nzdata"), use.names=FALSE)

    SparseArraySeed(ans_dim, ans_nzindex, ans_nzdata, ans_dimnames, check=FALSE)
}

setMethod("rbind", "SparseArraySeed",
    function(...)
    {
        abind_SparseArraySeed_objects(list(...), along=1L)
    }
)

setMethod("cbind", "SparseArraySeed",
    function(...)
    {
        abind_SparseArraySeed_objects(list(...), along=2L)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary isometric" array transformations
###
### A "unary isometric" array transformation is a transformation that returns
### an array-like object with the same dimensions as the input and where each
### element is the result of applying a function to the corresponding element
### in the input.
###
### Note that some "unary isometric" transformations preserve sparsity (e.g.
### is.na(), nchar(), round(), sqrt(), log1p(), etc...) and others don't
### (e.g. is.finite(), !, log(), etc..). We only implement the former.
###
### All the "unary isometric" array transformations implemented in this
### section return a SparseArraySeed object of the same dimensions as the
### input SparseArraySeed object.
###

.UNARY_ISO_OPS <- c("is.na", "is.infinite", "is.nan", "tolower", "toupper")

for (.Generic in .UNARY_ISO_OPS) {
    setMethod(.Generic, "SparseArraySeed",
        function(x)
        {
            GENERIC <- match.fun(.Generic)
            new_nzdata <- GENERIC(x@nzdata)
            BiocGenerics:::replaceSlots(x, nzdata=new_nzdata, check=FALSE)
        }
    )
}

setMethod("nchar", "SparseArraySeed",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
    {
        new_nzdata <- nchar(x@nzdata, type=type, allowNA=allowNA, keepNA=keepNA)
        BiocGenerics:::replaceSlots(x, nzdata=new_nzdata, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

setMethod("anyNA", "SparseArraySeed",
    function(x, recursive=FALSE) anyNA(x@nzdata, recursive=recursive)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which()
###

.nzindex_order <- function(nzindex)
    do.call(order, lapply(ncol(nzindex):1L, function(along) nzindex[ , along]))

setMethod("which", "SparseArraySeed",
    function(x, arr.ind=FALSE, useNames=TRUE)
    {
        if (!identical(useNames, TRUE))
            warning(wmsg("'useNames' is ignored when 'x' is ",
                         "a SparseArraySeed object or derivative"))
        if (!isTRUEorFALSE(arr.ind))
            stop(wmsg("'arr.ind' must be TRUE or FALSE"))
        idx1 <- which(x@nzdata)
        nzindex1 <- x@nzindex[idx1, , drop=FALSE]
        oo <- .nzindex_order(nzindex1)
        ans <- nzindex1[oo, , drop=FALSE]
        if (arr.ind)
            return(ans)
        Mindex2Lindex(ans, dim=dim(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max(), min(), range(), sum(), prod(), any(), all()
###

setMethod("Summary", "SparseArraySeed",
    function(x, ..., na.rm=FALSE)
    {
        GENERIC <- match.fun(.Generic)
        if (length(list(...)) != 0L)
            stop(wmsg(.Generic, "() method for SparseArraySeed objects ",
                      "only accepts a single object"))
        ## Whether 'x' contains zeros or not doesn't make a difference for
        ## sum() and any().
        if (.Generic %in% c("sum", "any"))
            return(GENERIC(x@nzdata, na.rm=na.rm))
        ## Of course a typical SparseArraySeed object "contains" zeros
        ## (i.e. it would contain zeros if we converted it to a dense
        ## representation with sparse2dense()). However, this is not
        ## guaranteed so we need to make sure to properly handle the case
        ## where it doesn't (admittedly unusual and definitely an inefficient
        ## way to represent dense data!)
        x_has_zeros <- length(x@nzdata) < length(x)
        if (!x_has_zeros)
            return(GENERIC(x@nzdata, na.rm=na.rm))
        x_type <- typeof(x@nzdata)
        if (.Generic == "all") {
            ## Mimic what 'all(sparse2dense(x))' would do.
            if (x_type == "double")
                warning("coercing argument of type 'double' to logical")
            return(FALSE)
        }
        zero <- vector(x_type, length=1L)
        GENERIC(zero, x@nzdata, na.rm=na.rm)
    }
)

### We override the "range" method defined above via the "Summary" method
### because we want to support the 'finite' argument like S3 method
### base::range.default() does.

### S3/S4 combo for range.SparseArraySeed
range.SparseArraySeed <- function(..., na.rm=FALSE, finite=FALSE)
{
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("range() method for SparseArraySeed objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    x_has_zeros <- length(x@nzdata) < length(x)
    if (!x_has_zeros)
        return(range(x@nzdata, na.rm=na.rm, finite=finite))
    zero <- vector(typeof(x@nzdata), length=1L)
    range(zero, x@nzdata, na.rm=na.rm, finite=finite)
}
### The signature of all the members of the S4 "Summary" group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that the S4 methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument before.
setMethod("range", "SparseArraySeed",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.SparseArraySeed(x, ..., na.rm=na.rm, finite=finite)

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

.mean_SparseArraySeed <- function(x, na.rm=FALSE)
{
    s <- sum(x@nzdata, na.rm=na.rm)
    nval <- length(x)
    if (na.rm)
        nval <- nval - sum(is.na(x@nzdata))
    s / nval
}

### S3/S4 combo for mean.SparseArraySeed
mean.SparseArraySeed <- function(x, na.rm=FALSE, ...)
    .mean_SparseArraySeed(x, na.rm=na.rm, ...)
setMethod("mean", "SparseArraySeed", .mean_SparseArraySeed)

