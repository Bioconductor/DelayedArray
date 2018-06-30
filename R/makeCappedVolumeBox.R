### =========================================================================
### makeCappedVolumeBox()
### -------------------------------------------------------------------------
###

### 'maxvol' is assumed to be a single integer > 0 and < 'prod(dim)'
.make_capped_volume_hypercube_box <- function(maxvol, dim)
{
    ans <- dim
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
    q <- get_RegularArrayGrid_dim(dim, ans + 1L) /
         get_RegularArrayGrid_dim(dim, ans)
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

### 'maxvol' is assumed to be a single integer > 0 and < 'prod(dim)'
.make_capped_volume_scale_box <- function(maxvol, dim)
{
    r <- (maxvol / prod(dim)) ^ (1 / length(dim))
    as.integer(r * dim)
}

### 'maxvol' is assumed to be a single integer > 0 and < 'prod(dim)'
.make_capped_volume_FDGF_box <- function(maxvol, dim)
{
    p <- cumprod(dim)
    w <- which(p <= maxvol)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N == 1L) {
        by <- maxvol
    } else {
        by <- maxvol %/% as.integer(p[[N - 1L]])
    }
    c(head(dim, n=N-1L), by, rep.int(1L, length(dim)-N))
}

.make_capped_volume_LDGF_box <- function(maxvol, dim)
{
    rev(.make_capped_volume_FDGF_box(maxvol, rev(dim)))
}

### Return the dimensions of a box that satisfies the following properties:
###   (1) Has a volume as close as possibe to (but not bigger than) 'maxvol'.
###   (2) Fits in the "constraining box" i.e. in the box of dimensions 'dim'.
###   (3) Has a shape that is as close as possible to the requested shape.
makeCappedVolumeBox <- function(maxvol, dim, shape=c("hypercube",
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

    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)

    shape <- match.arg(shape)

    if (maxvol == 0L)
        return(integer(length(dim)))

    if (maxvol >= prod(dim))
        return(dim)

    FUN <- switch(shape,
                  hypercube=.make_capped_volume_hypercube_box,
                  scale=.make_capped_volume_scale_box,
                  `first-dim-grows-first`=.make_capped_volume_FDGF_box,
                  `last-dim-grows-first`=.make_capped_volume_LDGF_box,
                  stop("unsupported 'shape'"))
    FUN(maxvol, dim)
}

### A capped-volume box related utility.
### If 'viewport_shape' is "first-dim-grows-first", return a linear grid.
makeRegularArrayGridOfCappedLengthViewports <-
    function(refdim, viewport_maxlen, viewport_shape=c("hypercube",
                                                       "scale",
                                                       "first-dim-grows-first",
                                                       "last-dim-grows-first"))
{
    spacings <- makeCappedVolumeBox(viewport_maxlen, refdim, viewport_shape)
    RegularArrayGrid(refdim, spacings)
}

