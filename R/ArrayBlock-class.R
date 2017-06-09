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

.validate_DIM_slot <- function(x, slot_name="DIM")
{
    DIM <- slot(x, slot_name)
    if (!is.integer(DIM))
        return(wmsg2(sprintf("'%s' slot must be an integer vector", slot_name)))
    if (length(DIM) == 0L)
        return(wmsg2(sprintf("'%s' slot cannot be empty", slot_name)))
    if (S4Vectors:::anyMissingOrOutside(DIM, 0L))
        return(wmsg2(sprintf("'%s' slot cannot contain negative or NA values",
                             slot_name)))
    TRUE
}

.validate_ranges_slot <- function(x)
{
    if (length(x@ranges) != length(x@DIM))
        return(wmsg2("'DIM' and 'ranges' slots must have the same length"))

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
        DIM_in1string <- paste0(object@DIM, collapse=" x ")
        cat(class(object), " object on a ", DIM_in1string, " array: ", sep="")
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BlockGrid objects
###

setClass("BlockGrid",
    contains="List",
    representation(
        DIM="integer",       # Dimensions of the array the grid is on.
        block_dim="integer"  # Dimensions of the first block in the grid.
    ),
    prototype(elementType="ArrayBlock")
)

.validate_BlockGrid <- function(x)
{
    msg <- .validate_DIM_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_DIM_slot(x, "block_dim")
    if (!isTRUE(msg))
        return(msg)
    if (length(x@DIM) != length(x@block_dim))
        return(wmsg2("'DIM' and 'block_dim' slots must have the same length"))
    if (!all(x@block_dim <= x@DIM))
        return(wmsg2("first block in the grid does not fit in the array ",
                     "the grid is on"))
    if (!any(x@DIM == 0L) && any(x@block_dim == 0L))
        return(wmsg2("first block in the grid cannot be empty unless",
                     "the grid is on an empty array"))
    TRUE
}

setValidity2("BlockGrid", .validate_BlockGrid)

### If 'block_dim' is omitted, returns a grid made of only block covering the
### whole array the grid is on.
BlockGrid <- function(dim, block_dim=dim)
    new("BlockGrid", DIM=dim, block_dim=block_dim)

### Return the number of blocks along each dimension.
.get_max_steps_along_each_dim <- function(x)
{
    mapply(function(D, d) {
               q <- D %/% d
               if (D %% d == 0L) q else q + 1L
           },
           x@DIM, x@block_dim)
}

setMethod("length", "BlockGrid",
    function(x)
    {
        if (any(x@DIM == 0L))
            return(0L)
        prod(.get_max_steps_along_each_dim(x))
    }
)

.to_array_index <- function(i, dim)
{
    ans <- integer(length(dim))
    for (along in seq_along(dim)) {
        d <- dim[[along]]
        ans[[along]] <- offset <- i %% d
        i <- (i - offset) %/% d
    }
    ans
}

### Return an ArrayBlock object.
setMethod("getListElement", "BlockGrid",
    function(x, i, exact=TRUE)
    {
        i <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                             error.if.nomatch=TRUE)
        max_steps_along_each_dim <- .get_max_steps_along_each_dim(x)
        block_offsets <- .to_array_index(i - 1L, max_steps_along_each_dim)
        block_offsets <- block_offsets * x@block_dim
        block_end <- pmin(x@DIM, block_offsets + x@block_dim)
        ArrayBlock(x@DIM, IRanges(block_offsets + 1L, block_end))
    }
)

setMethod("show", "BlockGrid",
    function(object)
    {
        DIM_in1string <- paste0(object@DIM, collapse=" x ")
        cat(class(object), " object with ", length(object), " blocks ",
            "on a ", DIM_in1string, " array:\n", sep="")
        for (i in seq_along(object)) {
            s <- make_string_from_ArrayBlock(object[[i]], with.brackets=TRUE)
            cat("[[", i, "]]: ", s, "\n", sep="")
        }
    }
)

