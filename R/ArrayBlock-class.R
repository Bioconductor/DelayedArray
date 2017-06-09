### =========================================================================
### ArrayBlock objects
### -------------------------------------------------------------------------


### We don't extend the IRanges class because we don't want to inherit the
### full Ranges API (most operations in that API would do the wrong thing on
### ArrayBlock objects).
setClass("ArrayBlock",
    representation(
        DIM="integer",    # Dimensions of the array this block belongs to.
        ranges="IRanges"  # Must be parallel to the 'DIM' slot.
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_DIM_slot <- function(x)
{
    if (!is.integer(x@DIM))
        return(wmsg2("'x@DIM' must be an integer vector"))
    x_ndim <- length(x@DIM)
    if (x_ndim == 0L)
        return(wmsg2("'x@DIM' cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x@DIM, 0L))
        return(wmsg2("'x@DIM' cannot contain negative or NA values"))
    TRUE
}

.validate_ranges_slot <- function(x)
{
    if (length(x@ranges) != length(x@DIM))
        return(wmsg2("'x@ranges' and 'x@DIM' must have the same length"))

    ## Check that the block is contained in the array it belongs to.
    x_start <- start(x@ranges)
    x_end <- end(x@ranges)
    if (!(all(x_start >= 1L) && all(x_end <= x@DIM)))
        return(wmsg2("object represents a block that is not ",
                     "contained in the array it belongs to"))

    ## A block cannot be empty.
    x_width <- width(x@ranges)
    if (any(x_width == 0L))
        return(wmsg2("a block cannot be empty"))
    TRUE
}

.validate_ArrayBlock <- function(x)
{
    msg <- .validate_DIM_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_ranges_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("ArrayBlock", .validate_ArrayBlock)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("dim", "ArrayBlock", function(x) width(x@ranges))

setMethod("length", "ArrayBlock", function(x) prod(dim(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### If 'ranges' is omitted, return a block that covers the whole array.
ArrayBlock <- function(dim, ranges=NULL)
{
    if (is.null(ranges))
        ranges <- IRanges(rep.int(1L, length(dim)), dim)
    new("ArrayBlock", DIM=dim, ranges=ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

make_string_from_ArrayBlock <- function(block, dimnames=NULL,
                                        with.brackets=FALSE)
{
    block_dim <- dim(block)
    ans <- as.character(block@ranges)
    ans[block_dim == block@DIM] <- ""
    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(block_dim) == length(dimnames))
        usename_idx <- which(block_dim == 1L &
                             block@DIM != 1L &
                             lengths(dimnames) != 0L)
        ans[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                        start(block@ranges)[usename_idx],
                                        SIMPLIFY=FALSE)
    }
    if (ans[[1L]] == "" && with.brackets)
        ans[[1L]] <- " "
    ans <- paste0(ans, collapse=", ")
    if (with.brackets)
        ans <- paste0("[", ans, "]")
    ans
}

setMethod("show", "ArrayBlock",
    function(object)
    {
        dim_in1string <- paste0(object@DIM, collapse=" x ")
        cat(class(object), " object on a ", dim_in1string, " array: ", sep="")
        s <- make_string_from_ArrayBlock(object, with.brackets=TRUE)
        cat(s, "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeNindexFromArrayBlock()
###

### Used in HDF5Array!
makeNindexFromArrayBlock <- function(block, expand.RangeNSBS=FALSE)
{
    block_dim <- dim(block)
    ndim <- length(block_dim)
    Nindex <- vector(mode="list", length=ndim)
    is_not_missing <- block_dim < block@DIM
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        block_starts <- start(block@ranges)
        block_ends <- end(block@ranges)
        is_width1 <- block_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
            function(i) {
                range_start <- block_starts[[i]]
                range_end <- block_ends[[i]]
                upper_bound <- block@DIM[[i]]
                new2("RangeNSBS", subscript=c(range_start, range_end),
                                  upper_bound=upper_bound,
                                  check=FALSE)
            }
        )
    }
    Nindex[expand_idx] <- as.list(block@ranges[expand_idx])
    Nindex
}

