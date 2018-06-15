### =========================================================================
### ArrayViewport and ArrayGrid objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayViewport objects
###

### We don't extend the IRanges class because we don't want to inherit the
### full Ranges API (most operations in that API would do the wrong thing on
### ArrayViewport objects).
setClass("ArrayViewport",
    contains="Array",
    representation(
        refdim="integer",  # Dimensions of "the reference array" i.e. the
                           # array on top of which the viewport is defined
                           # (a.k.a. "the underlying array").
        ranges="IRanges"   # Must be parallel to the 'refdim' slot.
    )
)

### Validity

.validate_ranges_slot <- function(x)
{
    x_ranges <- x@ranges
    x_refdim <- x@refdim
    if (length(x_ranges) != length(x_refdim))
        return(wmsg2("'ranges' and 'refdim' slots must have the same length"))

    ## Check that the viewport is contained in the reference array.
    x_start <- start(x_ranges)
    x_end <- end(x_ranges)
    if (!(all(x_start >= 1L) && all(x_end <= x_refdim)))
        return(wmsg2("object represents a viewport that is not ",
                     "within the bounds of the reference array"))

    ## A viewport cannot be longer than 2^31-1.
    x_dim <- width(x_ranges)
    if (prod(x_dim) > .Machine$integer.max)
        return(wmsg2("a viewport cannot be longer than .Machine$integer.max"))
    TRUE
}

.validate_ArrayViewport <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_ranges_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("ArrayViewport", .validate_ArrayViewport)

### Getters

setGeneric("refdim", function(x) standardGeneric("refdim"))

setMethod("refdim", "ArrayViewport", function(x) x@refdim)

setMethod("ranges", "ArrayViewport", function(x) x@ranges)

setMethod("start", "ArrayViewport", function(x) start(ranges(x)))

setMethod("width", "ArrayViewport", function(x) width(ranges(x)))

setMethod("end", "ArrayViewport", function(x) end(ranges(x)))

### 'width(x)' and 'dim(x)' are synonyms.
setMethod("dim", "ArrayViewport", function(x) width(ranges(x)))

### Constructor

### If 'ranges' is omitted, return a viewport that covers the whole
### reference array.
ArrayViewport <- function(refdim, ranges=NULL)
{
    if (is.null(ranges))
        ranges <- IRanges(rep.int(1L, length(refdim)), refdim)
    new("ArrayViewport", refdim=refdim, ranges=ranges)
}

### Show

make_string_from_ArrayViewport <- function(viewport, dimnames=NULL,
                                           as.2Dslice=FALSE,
                                           with.brackets=FALSE)
{
    if (!isTRUEorFALSE(as.2Dslice))
        stop("'as.2Dslice' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.brackets))
        stop("'with.brackets' must be TRUE or FALSE")

    viewport_ranges <- ranges(viewport)
    ans <- as.character(viewport_ranges)

    ## Place "blank" subscripts.
    viewport_dim <- dim(viewport)
    viewport_refdim <- refdim(viewport)
    useblank <- viewport_dim == viewport_refdim
    ndim <- length(viewport_dim)
    if (as.2Dslice && ndim >= 3L)
        useblank[3:ndim] <- FALSE
    ans[useblank] <- ""

    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(dimnames) == ndim)
        usename_idx <- which(!useblank &
                             viewport_dim == 1L &
                             lengths(dimnames) != 0L)
        ans[usename_idx] <- unlist(mapply(`[`,
                                          dimnames[usename_idx],
                                          start(viewport_ranges)[usename_idx],
                                          SIMPLIFY=FALSE,
                                          USE.NAMES=FALSE))
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
        refdim_in1string <- paste0(refdim(object), collapse=" x ")
        cat(class(object), " object on a ", refdim_in1string, " array: ",
            sep="")
        s <- make_string_from_ArrayViewport(object, with.brackets=TRUE)
        cat(s, "\n", sep="")
    }
)

### makeNindexFromArrayViewport()

### Used in HDF5Array!
makeNindexFromArrayViewport <- function(viewport, expand.RangeNSBS=FALSE)
{
    viewport_ranges <- ranges(viewport)
    viewport_dim <- dim(viewport)
    viewport_refdim <- refdim(viewport)
    ndim <- length(viewport_dim)
    Nindex <- vector(mode="list", length=ndim)
    is_not_missing <- viewport_dim < viewport_refdim
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        viewport_starts <- start(viewport_ranges)
        viewport_ends <- end(viewport_ranges)
        is_width1 <- viewport_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
            function(i) {
                range_start <- viewport_starts[[i]]
                range_end <- viewport_ends[[i]]
                upper_bound <- viewport_refdim[[i]]
                new2("RangeNSBS", subscript=c(range_start, range_end),
                                  upper_bound=upper_bound,
                                  check=FALSE)
            }
        )
    }
    Nindex[expand_idx] <- as.list(as(viewport_ranges[expand_idx],
                                     "CompressedIntegerList"))

    Nindex
}

### 2 utilities for extracting/replacing blocks from/in an array-like object

extract_block <- function(x, viewport)
{
    stopifnot(identical(dim(x), refdim(viewport)))
    Nindex <- makeNindexFromArrayViewport(viewport)
    subset_by_Nindex(x, Nindex)
}

### Return the modified array.
replace_block <- function(x, viewport, block)
{
    stopifnot(identical(dim(x), refdim(viewport)),
              identical(dim(viewport), dim(block)))
    Nindex <- makeNindexFromArrayViewport(viewport)
    replace_by_Nindex(x, Nindex, block)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayGrid objects
###
### An ArrayGrid object represents a grid on top of an array (called "the
### reference array" or "the underlying array"). The ArrayGrid class is a
### virtual class with 2 concrete subclasses, ArbitraryArrayGrid and
### RegularArrayGrid, for representing an arbitrarily-spaced or a
### regularly-spaced grid, respectively. The API we implement on these objects
### is divided into 3 groups of methods:
###   1) One special method:
###      - refdim(x): Return the dimensions of the reference array.
###   2) Methods from the array API (an ArrayGrid object can be seen as an
###      array of ArrayViewport objects that cover the reference array without
###      overlapping with each others):
###      - dim(x): Return the number of grid elements (i.e. viewports) along
###        each dimension of the reference array.
###      - x[[i_1, i_2, ..., i_n]]: Multi-dimensional double bracket
###        subsetting. Return an ArrayViewport object.
###   3) Methods from the list API:
###      - length(): Return the total number of grid elements.
###      - x[[i]]: Linear double bracket subsetting. Return an ArrayViewport
###        object.
### Groups 2) and 3) give these objects 2 semantics: array-like and list-like.
### Note that length() and "linear double bracket subsetting" are consistent
### with dim() and "multi-dimensional double bracket subsetting", respectively.
### So the array-like and list-like semantics are compatible.
###

setClass("ArrayGrid",
    contains=c("Array", "List"),
    representation("VIRTUAL"),
    prototype(elementType="ArrayViewport")
)

setClass("ArbitraryArrayGrid",
    contains="ArrayGrid",
    representation(
        tickmarks="list"      # A list of integer vectors, one along each
                              # dimension of the reference array,
                              # representing the tickmarks along that
                              # dimension. Each integer vector must be sorted
                              # in ascending order.
    )
)

setClass("RegularArrayGrid",
    contains="ArrayGrid",
    representation(
        refdim="integer",     # Dimensions of the reference array.
        spacings="integer"    # Grid spacing along each dimension.
    )
)

### Low-level helpers

.get_ArbitraryArrayGrid_spacings_along <- function(x, along)
    S4Vectors:::diffWithInitialZero(x@tickmarks[[along]])

.get_ArbitraryArrayGrid_max_spacings <- function(x)
{
    vapply(seq_along(x@tickmarks),
        function(along)
            max(0L, .get_ArbitraryArrayGrid_spacings_along(x, along)),
        integer(1),
        USE.NAMES=FALSE
    )
}

### Get length of biggest viewport in ArbitraryArrayGrid object 'x'.
.get_ArbitraryArrayGrid_maxlength <- function(x)
{
    ans <- prod(.get_ArbitraryArrayGrid_max_spacings(x))
    if (ans <= .Machine$integer.max)
        ans <- as.integer(ans)
    ans
}

### Get length of biggest viewport in RegularArrayGrid object 'x'.
.get_RegularArrayGrid_maxlength <- function(x)
{
    ans <- prod(x@spacings)
    if (ans <= .Machine$integer.max)
        ans <- as.integer(ans)
    ans
}

.get_RegularArrayGrid_dim <- function(refdim, spacings)
{
    ans <- refdim %/% spacings + (refdim %% spacings != 0L)
    ans[is.na(ans)] <- 1L
    ans
}

.get_RegularArrayGrid_spacings_along <- function(x, along)
{
    D <- x@refdim[[along]]
    if (D == 0L)
        return(0L)
    spacing <- x@spacings[[along]]
    ans <- rep.int(spacing, D %/% spacing)
    r <- D %% spacing
    if (r != 0L)
        ans <- c(ans, r)
    ans
}

### Validity

.valid_tickmarks <- function(tm)
{
    is.integer(tm) && !S4Vectors:::anyMissingOrOutside(tm, 0L) && isSorted(tm)
}
.validate_ArbitraryArrayGrid <- function(x)
{
    x_tickmarks <- x@tickmarks
    if (!is.list(x_tickmarks))
        return(wmsg2("'tickmarks' slot must be a list"))
    ok <- vapply(x_tickmarks, .valid_tickmarks, logical(1), USE.NAMES=FALSE)
    if (!all(ok))
        return(wmsg2("each list element in 'tickmarks' slot must be a ",
                     "sorted integer vector of non-negative values"))
    x_maxlen <- .get_ArbitraryArrayGrid_maxlength(x)
    if (x_maxlen > .Machine$integer.max)
        return(wmsg2("grid is too coarse (all grid elements must have a ",
                     "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("ArbitraryArrayGrid", .validate_ArbitraryArrayGrid)

.validate_RegularArrayGrid <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dim_slot(x, "spacings")
    if (!isTRUE(msg))
        return(msg)
    x_spacings <- x@spacings
    x_refdim <- x@refdim
    if (length(x_spacings) != length(x_refdim))
        return(wmsg2("'spacings' and 'refdim' slots must have ",
                     "the same length"))
    if (!all(x_spacings <= x_refdim))
        return(wmsg2("values in 'spacings' slot must be <= their ",
                     "corresponding value in 'refdim' slot"))
    if (any(x_spacings == 0L & x_refdim != 0L))
        return(wmsg2("values in 'spacings' slot cannot be 0 unless their ",
                     "corresponding value in 'refdim' slot is also 0"))
    x_maxlen <- .get_RegularArrayGrid_maxlength(x)
    if (x_maxlen > .Machine$integer.max)
        return(wmsg2("grid is too coarse (all grid elements must have a ",
                     "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("RegularArrayGrid", .validate_RegularArrayGrid)

### Getters

setMethod("refdim", "ArbitraryArrayGrid",
    function(x)
    {
        mapply(
            function(tm, tm_len) if (tm_len == 0L) 0L else tm[[tm_len]],
            x@tickmarks,
            lengths(x@tickmarks),
            USE.NAMES=FALSE
        )
    }
)

setMethod("refdim", "RegularArrayGrid", function(x) x@refdim)

setMethod("dim", "ArbitraryArrayGrid", function(x) lengths(x@tickmarks))

setMethod("dim", "RegularArrayGrid",
    function(x) .get_RegularArrayGrid_dim(refdim(x), x@spacings)
)

### Constructors

ArbitraryArrayGrid <- function(tickmarks)
{
    if (!is.list(tickmarks))
        stop(wmsg("'tickmarks' must be a list"))
    new("ArbitraryArrayGrid", tickmarks=tickmarks)
}

### Note that none of the dimensions of an RegularArrayGrid object can be 0,
### even when some dimensions of the reference array are 0 (in which case,
### the corresponding dimensions of the grid object are set to 1). As a
### consequence, an RegularArrayGrid object always contains at least 1 grid
### element. Each dimension of the first grid element is always equal to the
### spacing along that dimension i.e. for any RegularArrayGrid object,
### 'dim(grid[[1]])' is identical to 'spacings'.
### If 'spacings' is omitted, return a grid with a single grid element
### covering the whole reference array.
RegularArrayGrid <- function(refdim, spacings=refdim)
{
    if (!is.numeric(refdim))
        stop(wmsg("'refdim' must be an integer vector"))
    if (!is.integer(refdim))
        refdim <- as.integer(refdim)
    if (!is.numeric(spacings))
        stop(wmsg("'spacings' must be an integer vector"))
    if (!is.integer(spacings))
        spacings <- as.integer(spacings)
    if (length(refdim) != length(spacings))
        stop(wmsg("'refdim' and 'spacings' must have the same length"))
    new("RegularArrayGrid", refdim=refdim, spacings=spacings)
}

### [[

setMethod("getArrayElement", "ArbitraryArrayGrid",
    function(x, subscripts)
    {
        x_refdim <- refdim(x)
        ans_end <- mapply(`[[`, x@tickmarks, subscripts, USE.NAMES=FALSE)
        ans_width <- mapply(
            function(along, i)
                .get_ArbitraryArrayGrid_spacings_along(x, along)[[i]],
            seq_along(x_refdim),
            subscripts,
            USE.NAMES=FALSE
        )
        ans_ranges <- IRanges(end=ans_end, width=ans_width)
        ArrayViewport(x_refdim, ans_ranges)
    }
)

setMethod("getArrayElement", "RegularArrayGrid",
    function(x, subscripts)
    {
        x_refdim <- refdim(x)
        ans_offset <- (subscripts - 1L) * x@spacings
        ans_end <- pmin(ans_offset + x@spacings, refdim(x))
        ans_ranges <- IRanges(start=ans_offset + 1L, end=ans_end)
        ArrayViewport(x_refdim, ans_ranges)
    }
)

### dims() and lengths()

### NOT exported.
setGeneric("get_spacings_along", signature="x",
    function(x, along) standardGeneric("get_spacings_along")
)
setMethod("get_spacings_along", "ArbitraryArrayGrid",
    .get_ArbitraryArrayGrid_spacings_along
)
setMethod("get_spacings_along", "RegularArrayGrid",
    .get_RegularArrayGrid_spacings_along
)

### Equivalent to 't(vapply(x, dim, refdim(x)))' but faster.
setGeneric("dims", function(x) standardGeneric("dims"))

setMethod("dims", "ArrayGrid",
    function(x)
    {
        ans <- as.matrix(get_spacings_along(x, 1L))
        x_ndim <- length(refdim(x))
        if (x_ndim >= 2L) {
            for (along in 2:x_ndim) {
                spacings_along <- get_spacings_along(x, along)
                ans <- cbind(apply(ans, 2L, rep.int, length(spacings_along)),
                             rep(spacings_along, each=nrow(ans)))
            }
        }
        ans
    }
)

### Equivalent to 'vapply(x, length, integer(1))' or to 'rowProds(dims(x))'
### but faster.
### The sum of the hyper-volumes of all the grid elements should be equal
### to the hyper-volume of the reference array.
### More concisely: sum(lengths(x)) should be equal to 'prod(refdim(x))'.
setMethod("lengths", "ArrayGrid",
    function(x, use.names=TRUE)
    {
        ans <- get_spacings_along(x, 1L)
        x_ndim <- length(refdim(x))
        if (x_ndim >= 2L) {
            for (along in 2:x_ndim) {
                spacings_along <- get_spacings_along(x, along)
                ans <- ans * rep(spacings_along, each=length(ans))
            }
        }
        ans
    }
)

### Equivalent to 'max(lengths(x))' except that when 'x' is an ArrayGrid
### object this can be computed without computing 'lengths(x)' first so is
### very efficient.
setGeneric("maxlength", function(x) standardGeneric("maxlength"))
setMethod("maxlength", "ANY", function(x) max(lengths(x)))
setMethod("maxlength", "ArbitraryArrayGrid", .get_ArbitraryArrayGrid_maxlength)
setMethod("maxlength", "RegularArrayGrid", .get_RegularArrayGrid_maxlength)

### Show

### S3/S4 combo for as.character.ArrayGrid
.as.character.ArrayGrid <- function(x, with.brackets=FALSE)
{
    data <- vapply(x,
        function(viewport)
            make_string_from_ArrayViewport(viewport,
                                           with.brackets=with.brackets),
        character(1),
        USE.NAMES=FALSE
    )
    array(data, dim(x))
}
as.character.ArrayGrid <- function(x, ...) .as.character.ArrayGrid(x, ...)
setMethod("as.character", "ArrayGrid", .as.character.ArrayGrid)

setMethod("show", "ArrayGrid",
    function(object)
    {
        dim_in1string <- paste0(dim(object), collapse=" x ")
        refdim_in1string <- paste0(refdim(object), collapse=" x ")
        cat(dim_in1string, " ", class(object), " object ",
            "on a ", refdim_in1string, " array:\n", sep="")
        ## Turn 'object' into a character array.
        print(as.character(object, TRUE), quote=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGrid() generic and methods
###

### 'aind' must be a numeric vector or matrix (a vector is treated
### like a 1-row matrix).
setGeneric("mapToGrid", signature="grid",
    function(aind, grid, linear=FALSE) standardGeneric("mapToGrid")
)

.major_minor_as_list <- function(major, minor, grid, linear=FALSE)
{
    if (linear) {
        major <- linearInd(major, dim(grid))
        minor <- linearInd(minor, dims(grid)[major, , drop=FALSE])
    }
    list(major=major, minor=minor)
}

setMethod("mapToGrid", "ArbitraryArrayGrid",
    function(aind, grid, linear=FALSE)
    {
        if (!isTRUEorFALSE(linear))
            stop("'linear' must be TRUE or FALSE")
        ndim <- length(grid@tickmarks)
        aind <- normarg_aind(aind, ndim)
        major <- lapply(seq_len(ndim),
            function(along) {
                findInterval(aind[ , along], grid@tickmarks[[along]] + 1L) + 1L
            }
        )
        minor <- lapply(seq_len(ndim),
            function(along) {
                tm <- grid@tickmarks[[along]]
                tm_len <- length(tm)
                if (tm_len == 0L)
                    return(integer(0))
                offset <- c(0L, tm[-tm_len])
                aind[ , along] - offset[major[[along]]]
            }
        )
        major <- do.call(cbind, major)
        minor <- do.call(cbind, minor)
        .major_minor_as_list(major, minor, grid, linear=linear)
    }
)

setMethod("mapToGrid", "RegularArrayGrid",
    function(aind, grid, linear=FALSE)
    {
        if (!isTRUEorFALSE(linear))
            stop("'linear' must be TRUE or FALSE")
        ndim <- length(grid@spacings)
        aind <- normarg_aind(aind, ndim)
        d <- rep(grid@spacings, each=nrow(aind))
        aind0 <- aind - 1L  # 0-based indices
        major <- 1L + aind0 %/% d
        minor <- 1L + aind0 %% d
        .major_minor_as_list(major, minor, grid, linear=linear)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### downsample() generic and methods
###
### Reduce the "resolution" of a grid by the specified ratio.
### Act as an endomorphism.
###

setGeneric("downsample", signature="x",
    function(x, ratio=1L) standardGeneric("downsample")
)

.normarg_ratio <- function(ratio, x_dim)
{
    if (!is.numeric(ratio))
        stop(wmsg("'ratio' must be an integer vector"))
    if (!is.integer(ratio))
        ratio <- as.integer(ratio)
    ndim <- length(x_dim)
    if (length(ratio) != 1L && length(ratio) != ndim)
        stop(wmsg("'length(ratio)' must be 1 or the sane as 'length(dim(x))'"))
    if (S4Vectors:::anyMissingOrOutside(ratio, 0L))
        stop(wmsg("'ratio' cannot contain negative or NA values"))
    if (length(ratio) != ndim)
        ratio <- rep.int(ratio, ndim)
    if (any(ratio == 0L & x_dim != 0L))
        stop(wmsg("values in 'ratio' cannot be 0 unless their ",
                  "corresponding dimension in 'x' is also 0"))
    ratio
}

setMethod("downsample", "ArbitraryArrayGrid",
    function(x, ratio=1L)
    {
        ratio <- .normarg_ratio(ratio, dim(x))
        ans_tickmarks <- mapply(
            function(tm, tm_len, r) {
                if (tm_len == 0L)
                    return(integer(0))
                tm[seq2(tm_len, r)]
            },
            x@tickmarks,
            lengths(x@tickmarks),
            ratio,
            SIMPLIFY=FALSE,
            USE.NAMES=FALSE
        )
        ArbitraryArrayGrid(ans_tickmarks)
    }
)

setMethod("downsample", "RegularArrayGrid",
    function(x, ratio=1L)
    {
        ratio <- .normarg_ratio(ratio, dim(x))
        ans_spacings <- pmin(x@spacings * ratio, refdim(x))
        RegularArrayGrid(refdim(x), ans_spacings)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_spacings_for_hypercube_capped_length_blocks()
###
### Typically used to create a regular grid with a first block that is as
### close as possible to an hypercube and that has a length as close as
### possibe to (but not bigger than) 'block_maxlen'.
### NOT exported but used in HDF5Array!
get_spacings_for_hypercube_capped_length_blocks <-
    function(refdim, block_maxlen)
{
    if (!isSingleNumber(block_maxlen))
        stop("'block_maxlen' must be a single integer")
    if (!is.integer(block_maxlen))
        block_maxlen <- as.integer(block_maxlen)

    p <- prod(refdim)
    if (p <= block_maxlen)
        return(refdim)

    spacings <- refdim
    L <- max(spacings)
    while (TRUE) {
        is_max <- spacings == L
        not_max_spacings <- spacings[!is_max]
        L <- (block_maxlen / prod(not_max_spacings)) ^ (1 / sum(is_max))
        if (length(not_max_spacings) == 0L)
            break
        L2 <- max(not_max_spacings)
        if (L >= L2)
            break
        L <- L2
        spacings[is_max] <- L
    }
    spacings[is_max] <- as.integer(L)
    q <- .get_RegularArrayGrid_dim(refdim, spacings + 1L) /
         .get_RegularArrayGrid_dim(refdim, spacings)
    for (along in which(is_max)[order(q[is_max])]) {
        spacings[[along]] <- spacings[[along]] + 1L
        p <- prod(spacings)
        if (p == block_maxlen)
            break
        if (p > block_maxlen) {
            spacings[[along]] <- spacings[[along]] - 1L
            break
        }
    }
    spacings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Linear blocks
###
### An array block is "linear" if it would be made of elements contiguous in
### memory if the reference array was an ordinary R array (where the fastest
### changing dimension is the first one).
### A grid element is "linear" if it defines a linear block.
###

setGeneric("isLinear", function(x) standardGeneric("isLinear"))

setMethod("isLinear", "ArrayViewport",
    function(x)
    {
        x_width <- width(x)
        idx <- which(x_width != refdim(x))
        if (length(idx) == 0L)
            return(TRUE)
        all(tail(x_width, n=-idx[[1L]]) == 1L)
    }
)

### If the 1st grid element is linear, then they all are.
setMethod("isLinear", "ArrayGrid",
    function(x)
    {
        if (length(x) == 0L)
            return(TRUE)
        isLinear(x[[1L]])
    }
)

### Typically used to create a regular grid with linear blocks of length as
### close as possibe to (but not bigger than) 'block_maxlen'.
### NOT exported but used in HDF5Array!
get_spacings_for_linear_capped_length_blocks <- function(refdim, block_maxlen)
{
    if (!isSingleNumber(block_maxlen))
        stop("'block_maxlen' must be a single integer")
    if (!is.integer(block_maxlen))
        block_maxlen <- as.integer(block_maxlen)

    p <- cumprod(refdim)
    w <- which(p <= block_maxlen)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N > length(refdim))
        return(refdim)
    if (N == 1L) {
        by <- block_maxlen
    } else {
        by <- block_maxlen %/% as.integer(p[[N - 1L]])
    }
    c(head(refdim, n=N-1L), by, rep.int(1L, length(refdim)-N))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some unexported utilities
###

### Used in the DelayedMatrixStats package!
get_spacings_for_capped_length_blocks <-
    function(refdim, block_maxlen, block_shape=c("hypercube", "linear"))
{
    block_shape <- match.arg(block_shape)
    FUN <- switch(block_shape,
                  hypercube=get_spacings_for_hypercube_capped_length_blocks,
                  linear=get_spacings_for_linear_capped_length_blocks,
                  stop("unsupported 'block_shape'"))
    FUN(refdim, block_maxlen)
}

make_RegularArrayGrid_of_capped_length_blocks <-
    function(refdim, block_maxlen, block_shape=c("hypercube", "linear"))
{
    spacings <- get_spacings_for_capped_length_blocks(
                    refdim, block_maxlen, block_shape=block_shape)
    RegularArrayGrid(refdim, spacings)
}

### Used in unit tests.
split_array_in_capped_length_blocks <-
    function(x, block_maxlen, block_shape=c("hypercube", "linear"))
{
    grid <- make_RegularArrayGrid_of_capped_length_blocks(
                        dim(x), block_maxlen, block_shape=block_shape)
    lapply(grid, function(viewport) extract_block(x, viewport))
}

### Used in unit tests.
### Rebuild the original array from the blocks obtained with
### split_array_in_capped_length_blocks( , block_shape="hypercube") as
### an *ordinary* array. So if 'x' is an ordinary array, then:
###
###   blocks <- split_array_in_capped_length_blocks(x, block_maxlen,
###                                                 block_shape="hypercube")
###   unsplit_array_from_hypercube_blocks(blocks, x)
###
### should be a no-op for any 'block_maxlen' <= 'length(x)'.
unsplit_array_from_hypercube_blocks <- function(blocks, x)
{
    stop("Not ready yet")
}

### Used in unit tests.
### Rebuild the original array from the blocks obtained with
### split_array_in_capped_length_blocks( , block_shape="linear") as
### an *ordinary* array. So if 'x' is an ordinary array, then:
###
###   blocks <- split_array_in_capped_length_blocks(x, block_maxlen,
###                                                 block_shape="linear")
###   unsplit_array_from_linear_blocks(blocks, x)
###
### should be a no-op for any 'block_maxlen' <= 'length(x)'.
unsplit_array_from_linear_blocks <- function(blocks, x)
{
    ans <- combine_array_objects(blocks)
    set_dim(ans, dim(x))
}

