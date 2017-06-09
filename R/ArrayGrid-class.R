### =========================================================================
### ArrayViewport objects
### -------------------------------------------------------------------------


### We don't extend the IRanges class because we don't want to inherit the
### full Ranges API (most operations in that API would do the wrong thing on
### ArrayViewport objects).
setClass("ArrayViewport",
    representation(
        DIM="integer",    # Dimensions of the "reference array" i.e. the
                          # array the viewport belongs to.
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

    ## Check that the viewport is contained in the reference array.
    x_start <- start(x@ranges)
    x_end <- end(x@ranges)
    if (!(all(x_start >= 1L) && all(x_end <= x@DIM)))
        return(wmsg2("object represents a viewport that is not ",
                     "withing the bounds of the reference array"))

    ## A viewport cannot be empty.
    x_width <- width(x@ranges)
    if (any(x_width == 0L))
        return(wmsg2("a viewport cannot be empty"))
    TRUE
}

.validate_ArrayViewport <- function(x)
{
    msg <- .validate_DIM_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_ranges_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("ArrayViewport", .validate_ArrayViewport)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("dim", "ArrayViewport", function(x) width(x@ranges))

setMethod("length", "ArrayViewport", function(x) prod(dim(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### If 'ranges' is omitted, return a viewport that covers the whole
### reference array.
ArrayViewport <- function(dim, ranges=NULL)
{
    if (is.null(ranges))
        ranges <- IRanges(rep.int(1L, length(dim)), dim)
    new("ArrayViewport", DIM=dim, ranges=ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

make_string_from_ArrayViewport <- function(viewport, dimnames=NULL,
                                           with.brackets=FALSE)
{
    viewport_dim <- dim(viewport)
    ans <- as.character(viewport@ranges)
    ans[viewport_dim == viewport@DIM] <- ""
    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(viewport_dim) == length(dimnames))
        usename_idx <- which(viewport_dim == 1L &
                             viewport@DIM != 1L &
                             lengths(dimnames) != 0L)
        ans[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                        start(viewport@ranges)[usename_idx],
                                        SIMPLIFY=FALSE)
    }
    if (ans[[1L]] == "" && with.brackets)
        ans[[1L]] <- " "
    ans <- paste0(ans, collapse=", ")
    if (with.brackets)
        ans <- paste0("[", ans, "]")
    ans
}

setMethod("show", "ArrayViewport",
    function(object)
    {
        DIM_in1string <- paste0(object@DIM, collapse=" x ")
        cat(class(object), " object on a ", DIM_in1string, " array: ", sep="")
        s <- make_string_from_ArrayViewport(object, with.brackets=TRUE)
        cat(s, "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeNindexFromArrayViewport()
###

### Used in HDF5Array!
makeNindexFromArrayViewport <- function(viewport, expand.RangeNSBS=FALSE)
{
    viewport_dim <- dim(viewport)
    ndim <- length(viewport_dim)
    Nindex <- vector(mode="list", length=ndim)
    is_not_missing <- viewport_dim < viewport@DIM
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        viewport_starts <- start(viewport@ranges)
        viewport_ends <- end(viewport@ranges)
        is_width1 <- viewport_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
            function(i) {
                range_start <- viewport_starts[[i]]
                range_end <- viewport_ends[[i]]
                upper_bound <- viewport@DIM[[i]]
                new2("RangeNSBS", subscript=c(range_start, range_end),
                                  upper_bound=upper_bound,
                                  check=FALSE)
            }
        )
    }
    Nindex[expand_idx] <- as.list(viewport@ranges[expand_idx])
    Nindex
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayGrid objects
###

setClass("ArrayGrid",
    contains="List",
    representation(
        DIM="integer",       # Dimensions of the array the grid is on.
        viewport_dim="integer"  # Dimensions of the first viewport in the grid.
    ),
    prototype(elementType="ArrayViewport")
)

.validate_ArrayGrid <- function(x)
{
    msg <- .validate_DIM_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_DIM_slot(x, "viewport_dim")
    if (!isTRUE(msg))
        return(msg)
    if (length(x@DIM) != length(x@viewport_dim))
        return(wmsg2("'DIM' and 'viewport_dim' slots must have ",
                     "the same length"))
    if (!all(x@viewport_dim <= x@DIM))
        return(wmsg2("first viewport in the grid does not fit in the ",
                     "reference array"))
    if (any(x@DIM != 0L & x@viewport_dim == 0L))
        return(wmsg2("first viewport in the grid cannot have any of its ",
                     "dimensions set to 0 unless the corresponding ",
                     "dimension in the reference array is set to 0"))
    TRUE
}

setValidity2("ArrayGrid", .validate_ArrayGrid)

### If 'viewport_dim' is omitted, return a grid made of a single viewport
### covering the whole reference array.
ArrayGrid <- function(dim, viewport_dim=dim)
    new("ArrayGrid", DIM=dim, viewport_dim=viewport_dim)

### Return the number of viewports along each dimension of the reference
### array.
setMethod("dim", "ArrayGrid",
    function(x)
    {
        mapply(function(D, d) {
                   if (d <= 0L) return(0L)
                   q <- D %/% d
                   if (D %% d == 0L) q else q + 1L
               },
               x@DIM,
               x@viewport_dim)
    }
)

setMethod("length", "ArrayGrid", function(x) prod(dim(x)))

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

### Return an ArrayViewport object.
setMethod("getListElement", "ArrayGrid",
    function(x, i, exact=TRUE)
    {
        i <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                             error.if.nomatch=TRUE)
        viewport_offsets <- .to_array_index(i - 1L, dim(x))
        viewport_offsets <- viewport_offsets * x@viewport_dim
        viewport_ends <- pmin(x@DIM, viewport_offsets + x@viewport_dim)
        ArrayViewport(x@DIM, IRanges(viewport_offsets + 1L, viewport_ends))
    }
)

setMethod("show", "ArrayGrid",
    function(object)
    {
        DIM_in1string <- paste0(object@DIM, collapse=" x ")
        cat(class(object), " object with ", length(object), " viewports ",
            "on a ", DIM_in1string, " array:\n", sep="")
        for (i in seq_along(object)) {
            s <- make_string_from_ArrayViewport(object[[i]], with.brackets=TRUE)
            cat("[[", i, "]]: ", s, "\n", sep="")
        }
    }
)

