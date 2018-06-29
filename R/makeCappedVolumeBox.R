### =========================================================================
### makeCappedVolumeBox()
### -------------------------------------------------------------------------
###


.make_capped_volume_hypercube_box <- function(maxvol, dim)
{
    refvol <- prod(dim)
    if (refvol <= maxvol)
        return(dim)

    spacings <- dim
    L <- max(spacings)
    while (TRUE) {
        is_max <- spacings == L
        not_max_spacings <- spacings[!is_max]
        L <- (maxvol / prod(not_max_spacings)) ^ (1 / sum(is_max))
        if (length(not_max_spacings) == 0L)
            break
        L2 <- max(not_max_spacings)
        if (L >= L2)
            break
        L <- L2
        spacings[is_max] <- L
    }
    spacings[is_max] <- as.integer(L)
    q <- get_RegularArrayGrid_dim(dim, spacings + 1L) /
         get_RegularArrayGrid_dim(dim, spacings)
    for (along in which(is_max)[order(q[is_max])]) {
        spacings[[along]] <- spacings[[along]] + 1L
        p <- prod(spacings)
        if (p == maxvol)
            break
        if (p > maxvol) {
            spacings[[along]] <- spacings[[along]] - 1L
            break
        }
    }
    spacings
}

.make_capped_volume_proportional_box <- function(maxvol, dim)
{
    stop("'shape=\"proportional\"'")
}

.make_capped_volume_linear_box <- function(maxvol, dim)
{
    p <- cumprod(dim)
    w <- which(p <= maxvol)
    N <- if (length(w) == 0L) 1L else w[[length(w)]] + 1L
    if (N > length(dim))
        return(dim)
    if (N == 1L) {
        by <- maxvol
    } else {
        by <- maxvol %/% as.integer(p[[N - 1L]])
    }
    c(head(dim, n=N-1L), by, rep.int(1L, length(dim)-N))
}

### Return the dimensions of a box that satisfies the following properties:
###   (1) Has a volume as close as possibe to (but not bigger than) 'maxvol'.
###   (2) Fits in the "reference box" i.e. in the box of dimensions 'dim'.
###   (3) Has a shape that is as close as possible to the requested shape.
makeCappedVolumeBox <- function(maxvol, dim,
                                shape=c("hypercube", "proportional", "linear"))
{
    if (!isSingleNumber(maxvol))
        stop("'maxvol' must be a single integer")
    if (!is.integer(maxvol))
        maxvol <- as.integer(maxvol)
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    shape <- match.arg(shape)
    FUN <- switch(shape,
                  hypercube=.make_capped_volume_hypercube_box,
                  proportional=.make_capped_volume_proportional_box,
                  linear=.make_capped_volume_linear_box,
                  stop("unsupported 'shape'"))
    FUN(maxvol, dim)
}


### A capped-volume box related utility.
makeRegularArrayGridOfCappedLengthViewports <-
    function(refdim, viewport_maxlen,
             viewport_shape=c("hypercube", "proportional", "linear"))
{
    spacings <- makeCappedVolumeBox(viewport_maxlen, refdim, viewport_shape)
    RegularArrayGrid(refdim, spacings)
}

