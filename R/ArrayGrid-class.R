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
    representation(
        refdim="integer",  # Dimensions of "the reference array" i.e. the
                           # array on top of which the viewport is defined
                           # (a.k.a. "the underlying array").
        ranges="IRanges"   # Must be parallel to the 'refdim' slot.
    )
)

### Validity

.validate_refdim_slot <- function(x, slotname="refdim")
{
    x_refdim <- slot(x, slotname)
    if (!is.integer(x_refdim))
        return(wmsg2(sprintf("'%s' slot must be an integer vector", slotname)))
    if (length(x_refdim) == 0L)
        return(wmsg2(sprintf("'%s' slot cannot be empty", slotname)))
    if (S4Vectors:::anyMissingOrOutside(x_refdim, 0L))
        return(wmsg2(sprintf("'%s' slot cannot contain negative or NA values",
                             slotname)))
    TRUE
}

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
                     "withing the bounds of the reference array"))

    ## A viewport cannot be empty.
    x_width <- width(x_ranges)
    if (any(x_width == 0L))
        return(wmsg2("a viewport cannot be empty"))

    ## A viewport cannot be longer than 2^31-1.
    if (prod(x_width) > .Machine$integer.max)
        return(wmsg2("a vewport cannot be longer than .Machine$integer.max"))
    TRUE
}

.validate_ArrayViewport <- function(x)
{
    msg <- .validate_refdim_slot(x)
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

### This is the hyper-volume of the viewport.
setMethod("length", "ArrayViewport", function(x) prod(dim(x)))

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
                                           with.brackets=FALSE)
{
    if (!isTRUEorFALSE(with.brackets))
        stop("'with.brackets' must be TRUE or FALSE")
    viewport_ranges <- ranges(viewport)
    viewport_dim <- dim(viewport)
    viewport_refdim <- refdim(viewport)
    ans <- as.character(viewport_ranges)
    ans[viewport_dim == viewport_refdim] <- ""
    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(viewport_dim) == length(dimnames))
        usename_idx <- which(viewport_dim == 1L &
                             viewport_refdim != 1L &
                             lengths(dimnames) != 0L)
        ans[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                        start(viewport_ranges)[usename_idx],
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
    Nindex[expand_idx] <- as.list(viewport_ranges[expand_idx])
    Nindex
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayGrid objects
###
### An ArrayGrid object represents a grid on top of an array (called "the
### reference array" or "the underlying array"). The ArrayGrid class is a
### virtual class with 2 concrete subclasses, ArrayArbitraryGrid and
### ArrayRegularGrid, for representing an arbitrarily-spaced or a
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
    contains="List",
    representation("VIRTUAL"),
    prototype(elementType="ArrayViewport")
)

setClass("ArrayArbitraryGrid",
    contains="ArrayGrid",
    representation(
        tickmarks="list"      # A list of integer vectors, one along each
                              # dimension of the reference array,
                              # representing the tickmarks along that
                              # dimension. Each integer vector must be sorted
                              # in strictly ascending order.
    )
)

setClass("ArrayRegularGrid",
    contains="ArrayGrid",
    representation(
        refdim="integer",     # Dimensions of the reference array.
        spacings="integer"    # Dimensions of the first viewport in the grid.
    )
)

### Low-level helpers

.get_ArrayArbitraryGrid_spacings_along <- function(x, along)
    S4Vectors:::diffWithInitialZero(x@tickmarks[[along]])

.get_ArrayArbitraryGrid_max_spacings <- function(x)
{
    vapply(seq_along(x@tickmarks),
           function(along)
               max(0L, .get_ArrayArbitraryGrid_spacings_along(x, along)),
           integer(1))
}

.get_ArrayRegularGrid_dim <- function(refdim, spacings)
{
    ans <- refdim %/% spacings + (refdim %% spacings != 0L)
    ans[is.na(ans)] <- 0L
    ans
}

.get_ArrayRegularGrid_spacings_along <- function(x, along)
{
    D <- x@refdim[[along]]
    spacing <- x@spacings[[along]]
    ans <- rep.int(spacing, D %/% spacing)
    r <- D %% spacing
    if (r != 0L)
        ans <- c(ans, r)
    ans
}

### Validity

.validate_tickmarks <- function(tm)
{
    is.integer(tm) &&
        !S4Vectors:::anyMissingOrOutside(tm, 1L) &&
        isStrictlySorted(tm)
}
.validate_ArrayArbitraryGrid <- function(x)
{
    x_tickmarks <- x@tickmarks
    if (!is.list(x_tickmarks))
        return(wmsg2("'tickmarks' slot must be a list"))
    ok <- vapply(x_tickmarks, .validate_tickmarks, logical(1), USE.NAMES=FALSE)
    if (!all(ok))
        return(wmsg2("each list element in 'tickmarks' slot must be an ",
                     "integer vector of strictly ascending positive values"))
    max_spacings <- .get_ArrayArbitraryGrid_max_spacings(x)
    if (prod(max_spacings) > .Machine$integer.max)
        return(wmsg2("grid is too coarse (all grid elements must have a ",
                     "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("ArrayArbitraryGrid", .validate_ArrayArbitraryGrid)

.validate_ArrayRegularGrid <- function(x)
{
    msg <- .validate_refdim_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_refdim_slot(x, "spacings")
    if (!isTRUE(msg))
        return(msg)
    x_spacings <- x@spacings
    x_refdim <- x@refdim
    if (length(x_spacings) != length(x_refdim))
        return(wmsg2("'spacings' and 'refdim' slots must have ",
                     "the same length"))
    if (!all(x_spacings <= x_refdim))
        return(wmsg2("first viewport in the grid does not fit in the ",
                     "reference array"))
    if (any(x_spacings == 0L & x_refdim != 0L))
        return(wmsg2("first viewport in the grid cannot have any of its ",
                     "dimensions set to 0 unless the corresponding ",
                     "dimension in the reference array is set to 0"))
    if (prod(x_spacings) > .Machine$integer.max)
        return(wmsg2("grid is too coarse (all grid elements must have a ",
                     "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("ArrayRegularGrid", .validate_ArrayRegularGrid)

### Getters

setMethod("refdim", "ArrayArbitraryGrid",
    function(x)
    {
        mapply(function(tm, tm_len) if (tm_len == 0L) 0L else tm[[tm_len]],
               x@tickmarks,
               lengths(x@tickmarks),
               USE.NAMES=FALSE)
    }
)

setMethod("refdim", "ArrayRegularGrid", function(x) x@refdim)

setMethod("dim", "ArrayArbitraryGrid", function(x) lengths(x@tickmarks))

setMethod("dim", "ArrayRegularGrid",
    function(x) .get_ArrayRegularGrid_dim(refdim(x), x@spacings)
)

setMethod("length", "ArrayGrid", function(x) prod(dim(x)))

### Constructors

ArrayArbitraryGrid <- function(tickmarks)
    new("ArrayArbitraryGrid", tickmarks=tickmarks)

### If 'spacings' is omitted, return a grid made of a single viewport covering
### the whole reference array.
ArrayRegularGrid <- function(refdim, spacings=refdim)
    new("ArrayRegularGrid", refdim=refdim, spacings=spacings)

### [[

### Implement multi-dimensional double bracket subsetting.
### 'subscripts' is assumed to be an integer vector parallel to 'dim(x)' and
### with no out-of-bounds subscripts (i.e. 'all(subscripts >= 1)' and
### 'all(subscripts <= dim(x))').
### NOT exported for now but should probably be at some point (like
### S4Vectors::getListElement() is).
setGeneric("getArrayElement", signature="x",
    function(x, subscripts) standardGeneric("getArrayElement")
)

setMethod("getArrayElement", "ArrayArbitraryGrid",
    function(x, subscripts)
    {
        x_refdim <- refdim(x)
        ans_end <- mapply(`[[`, x@tickmarks, subscripts)
        ans_width <- mapply(
            function(along, i)
                .get_ArrayArbitraryGrid_spacings_along(x, along)[[i]],
            seq_along(x_refdim),
            subscripts)
        ans_ranges <- IRanges(end=ans_end, width=ans_width)
        ArrayViewport(x_refdim, ans_ranges)
    }
)

setMethod("getArrayElement", "ArrayRegularGrid",
    function(x, subscripts)
    {
        x_refdim <- refdim(x)
        ans_offset <- (subscripts - 1L) * x@spacings
        ans_end <- pmin(ans_offset + x@spacings, refdim(x))
        ans_ranges <- IRanges(start=ans_offset + 1L, end=ans_end)
        ArrayViewport(x_refdim, ans_ranges)
    }
)

### Return an integer vector parallel to 'dim' and guaranteed to contain no
### out-of-bounds subscripts.
.from_linear_to_multi_subscript <- function(i, dim)
{
    stopifnot(isSingleInteger(i))
    if (i < 1L || i > prod(dim))
        stop("subscript is out of bounds")
    i <- i - 1L
    subscripts <- integer(length(dim))
    for (along in seq_along(dim)) {
        d <- dim[[along]]
        subscripts[[along]] <- offset <- i %% d
        i <- (i - offset) %/% d
    }
    subscripts + 1L
}

### Support multi-dimensional and linear subsetting.
setMethod("[[", "ArrayGrid",
    function(x, i, j, ...)
    {
        if (missing(x))
            stop("'x' is missing")
        Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
        nsubscript <- length(Nindex)
        x_dim <- dim(x)
        x_ndim <- length(x_dim)
        if (!(nsubscript == 1L || nsubscript == x_ndim))
            stop("incorrect number of subscripts")
        ok <- vapply(Nindex, isSingleInteger, logical(1), USE.NAMES=FALSE)
        if (!all(ok))
            stop(wmsg("each subscript must be a single integer ",
                      "when subsetting an ArrayGrid object with [["))
        subscripts <- unlist(Nindex, use.names=FALSE)
        if (nsubscript != x_ndim) {
            ## Translate linear subsetting into multi-dimensional subsetting.
            subscripts <- .from_linear_to_multi_subscript(subscripts, x_dim)
        } else if (!(all(subscripts >= 1L) && all(subscripts <= x_dim))) {
            stop("some subscripts are out of bounds")
        }
        getArrayElement(x, subscripts)
    }
)

### lengths()

### NOT exported.
setGeneric("get_spacings_along", signature="x",
    function(x, along) standardGeneric("get_spacings_along")
)
setMethod("get_spacings_along", "ArrayArbitraryGrid",
    .get_ArrayArbitraryGrid_spacings_along
)
setMethod("get_spacings_along", "ArrayRegularGrid",
    .get_ArrayRegularGrid_spacings_along
)

### Equivalent to 'vapply(x, length, integer(1))' but faster.
### The sum of the hyper-volumes of all the grid elements should be equal
### to the hyper-volume of the reference array.
### More concisely: sum(lengths(x)) should be equal to 'prod(refdim(x))'.
setMethod("lengths", "ArrayGrid",
    function (x, use.names=TRUE)
    {
        ans <- get_spacings_along(x, 1L)
        x_ndim <- length(refdim(x))
        if (x_ndim >= 2L) {
            for (along in 2:x_ndim)
              ans <- ans * rep(get_spacings_along(x, along), each=length(ans))
        }
        ans
    }
)

### Show

### S3/S4 combo for as.character.ArrayGrid
.as.character.ArrayGrid <- function(x, with.brackets=FALSE)
{
    data <- vapply(x,
        function(viewport)
            make_string_from_ArrayViewport(viewport,
                                           with.brackets=with.brackets),
        character(1)
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
### A convenience helper.
###

extract_array_block <- function(x, grid, b)
{
    Nindex <- makeNindexFromArrayViewport(grid[[b]])
    subset_by_Nindex(x, Nindex)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_max_spacings_for_hypercube_blocks()
###

### Typically used to create a grid that divides the reference array 'x' into
### blocks that have a shape that is as close as possible to an hypercube
### (while having their length <= 'max_block_len'):
###
###    max_block_len <- get_max_block_length(type(x))
###    spacings <- get_max_spacings_for_hypercube_blocks(dim(x), max_block_len)
###    grid <- ArrayRegularGrid(dim(x), spacings)
###
### NOT exported but used in HDF5Array!
get_max_spacings_for_hypercube_blocks <- function(refdim, max_block_len)
{
    if (!isSingleNumber(max_block_len))
        stop("'max_block_len' must be a single number")
    p <- prod(refdim)
    if (p <= max_block_len)
        return(refdim)

    spacings <- refdim
    L <- max(spacings)
    while (TRUE) {
        is_max <- spacings == L
        not_max_spacings <- spacings[!is_max]
        L <- (max_block_len / prod(not_max_spacings)) ^ (1 / sum(is_max))
        if (length(not_max_spacings) == 0L)
            break
        L2 <- max(not_max_spacings)
        if (L >= L2)
            break
        L <- L2
        spacings[is_max] <- L
    }
    spacings[is_max] <- as.integer(L)
    q <- .get_ArrayRegularGrid_dim(refdim, spacings + 1L) /
         .get_ArrayRegularGrid_dim(refdim, spacings)
    for (along in which(is_max)[order(q[is_max])]) {
        spacings[[along]] <- spacings[[along]] + 1L
        p <- prod(spacings)
        if (p == max_block_len)
            break
        if (p > max_block_len) {
            spacings[[along]] <- spacings[[along]] - 1L
            break
        }
    }
    spacings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Grids with "linear grid elements"
###
### The grid elements are said to be "linear" if they divide the reference
### array in "linear blocks", i.e. in blocks that would be made of array
### elements contiguous in memory if the reference array was an ordinary R
### array (where the fastest changing dimension is the first one).
### Note that if the 1st grid element is linear, then they all are.
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

setMethod("isLinear", "ArrayGrid",
    function(x)
    {
        if (length(x) == 0L)
            return(TRUE)
        isLinear(x[[1L]])
    }
)

### Typically used to create a regular grid with "linear grid elements" i.e.
### a grid that divides the reference array 'x' into "linear blocks":
###
###    max_block_len <- get_max_block_length(type(x))
###    spacings <- get_max_spacings_for_linear_blocks(dim(x), max_block_len)
###    grid <- ArrayRegularGrid(dim(x), spacings)
###
### All the grid elements are guaranteed to have a length <= 'max_block_len'.
### NOT exported but used in HDF5Array!
get_max_spacings_for_linear_blocks <- function(refdim, max_block_len)
{
    if (!isSingleNumber(max_block_len))
        stop("'max_block_len' must be a single number")
    if (!is.integer(max_block_len))
        max_block_len <- as.integer(max_block_len)
    p <- cumprod(refdim)
    w <- which(p <= max_block_len)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N > length(refdim))
        return(refdim)
    if (N == 1L) {
        by <- max_block_len
    } else {
        by <- max_block_len %/% as.integer(p[[N - 1L]])
    }
    c(head(refdim, n=N-1L), by, rep.int(1L, length(refdim)-N))
}

### NOT exported but used in unit tests.
split_array_in_linear_blocks <- function(x, max_block_len)
{
    spacings <- get_max_spacings_for_linear_blocks(dim(x), max_block_len)
    grid <- ArrayRegularGrid(dim(x), spacings)
    lapply(seq_along(grid),
           function(b) extract_array_block(x, grid, b))
}

### NOT exported but used in unit tests.
### Rebuild the original array from the blocks obtained by
### split_array_in_linear_blocks() as an *ordinary* array.
### So if 'x' is an ordinary array, then:
###
###   blocks <- split_array_in_linear_blocks(x, max_block_len)
###   unsplit_array_from_linear_blocks(blocks, x)
###
### should be a no-op for any 'max_block_len' < 'length(x)'.
unsplit_array_from_linear_blocks <- function(blocks, x)
{
    ans <- combine_array_objects(blocks)
    dim(ans) <- dim(x)
    ans
}

