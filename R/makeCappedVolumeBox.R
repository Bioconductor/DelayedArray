### =========================================================================
### Utilities to make capped volume boxes
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeCappedVolumeBox()
###

### 'maxvol' is assumed to be a single integer >= 2 and < 'prod(maxdim)'.
.make_capped_volume_hypercube_box <- function(maxvol, maxdim)
{
    ans <- maxdim
    L <- max(ans)
    while (TRUE) {
        is_max <- ans == L
        not_max_ans <- ans[!is_max]
        L <- (maxvol / prod(not_max_ans)) ^ (1 / sum(is_max))
        if (length(not_max_ans) == 0L)
            break
        L2 <- max(not_max_ans)
        if (L >= L2)
            break
        L <- L2
        ans[is_max] <- L
    }
    ans[is_max] <- as.integer(L)
    q <- get_RegularArrayGrid_dim(maxdim, ans + 1L) /
         get_RegularArrayGrid_dim(maxdim, ans)
    for (along in which(is_max)[order(q[is_max])]) {
        ans[[along]] <- ans[[along]] + 1L
        p <- prod(ans)
        if (p == maxvol)
            break
        if (p > maxvol) {
            ans[[along]] <- ans[[along]] - 1L
            break
        }
    }
    ans
}

### 'maxvol' is assumed to be a single integer >= 2 and < 'prod(maxdim)'.
### The algo used below could be improved. For exampe it does some weird
### things like:
###     > .make_capped_volume_scale_box(11, c(3, 50, 10))
###     [1] 1 9 1
###     > .make_capped_volume_scale_box(12, c(3, 50, 10))
###     [1] 1 8 1
.make_capped_volume_scale_box <- function(maxvol, maxdim)
{
    ## Some good properties of shrinkbox():
    ## - The output dims are always >= 1.
    ## - If r is < 1, then input dims that are > 1 will decrease and those
    ##   at 1 will remain at 1.
    shrinkbox <- function(bdim, r) pmax(as.integer(bdim * r), 1L)

    p <- 1 / length(maxdim)
    bdim <- maxdim                  # all(maxdim >= 1) is TRUE
    ## Loop will typically go thru 2 to 18 iterations before it breaks.
    ## An example that requires 18 iterations:
    ## - maxvol <- 70000
    ## - maxdim <- c(30, 15000000)
    while (TRUE) {
        bvol <- prod(bdim)          # can't be 0
        if (bvol <= maxvol)
            break
        r <- (maxvol / bvol)^p      # < 1
        bdim <- shrinkbox(bdim, r)  # reduce all dims (except those already
                                    # at 1) so volume is guaranteed to reduce
                                    # at each loop
    }
    bdim
}

### 'maxvol' is assumed to be a single integer >= 2 and < 'prod(maxdim)'.
.make_capped_volume_FDGF_box <- function(maxvol, maxdim)
{
    p <- cumprod(maxdim)
    w <- which(p <= maxvol)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N == 1L) {
        by <- maxvol
    } else {
        by <- maxvol %/% as.integer(p[[N - 1L]])
    }
    c(head(maxdim, n=N-1L), by, rep.int(1L, length(maxdim)-N))
}

.make_capped_volume_LDGF_box <- function(maxvol, maxdim)
{
    rev(.make_capped_volume_FDGF_box(maxvol, rev(maxdim)))
}

### Return the dimensions of a box that satisfies the following properties:
###   1. Has a volume as close as possibe to (but not bigger than) 'maxvol'.
###   2. Fits in the "constraining box" i.e. in the box of dimensions 'maxdim'.
###   3. Has a non-zero volume if the "constraining box" has a non-zero volume.
###   4. Has a shape that is as close as possible to the requested shape.
makeCappedVolumeBox <- function(maxvol, maxdim, shape=c("hypercube",
                                                        "scale",
                                                        "first-dim-grows-first",
                                                        "last-dim-grows-first"))
{
    if (!isSingleNumber(maxvol))
        stop("'maxvol' must be a single integer")
    if (!is.integer(maxvol))
        maxvol <- as.integer(maxvol)
    if (maxvol < 0L)
        stop("'maxvol' must be a non-negative integer")

    if (!is.numeric(maxdim))
        stop(wmsg("'maxdim' must be an integer vector"))
    if (!is.integer(maxdim))
        maxdim <- as.integer(maxdim)

    shape <- match.arg(shape)

    if (maxvol == 0L || any(maxdim == 0L))
        return(integer(length(maxdim)))

    if (maxvol == 1L)
        return(rep.int(1L, length(maxdim)))

    if (maxvol >= prod(maxdim))
        return(maxdim)

    FUN <- switch(shape,
                  hypercube=.make_capped_volume_hypercube_box,
                  scale=.make_capped_volume_scale_box,
                  `first-dim-grows-first`=.make_capped_volume_FDGF_box,
                  `last-dim-grows-first`=.make_capped_volume_LDGF_box,
                  stop("unsupported 'shape'"))
    FUN(maxvol, maxdim)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeRegularArrayGridOfCappedLengthViewports()
###

### A capped-volume box related utility.
### If 'viewport_shape' is "first-dim-grows-first", return a linear grid.
makeRegularArrayGridOfCappedLengthViewports <-
    function(refdim, viewport_len, viewport_shape=c("hypercube",
                                                       "scale",
                                                       "first-dim-grows-first",
                                                       "last-dim-grows-first"))
{
    spacings <- makeCappedVolumeBox(viewport_len, refdim, viewport_shape)
    RegularArrayGrid(refdim, spacings)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Linear viewports and grids
###
### An array viewport is "linear" if it is made of reference array elements
### that would be contiguous in memory if the reference array was an ordinary
### R array (where the fastest changing dimension is the first one).
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

