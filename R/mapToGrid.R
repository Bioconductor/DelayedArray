### =========================================================================
### Map reference array positions to grid positions and vice-versa
### -------------------------------------------------------------------------


setGeneric("mapToGrid", signature="grid",
    function(Mindex, grid, linear=FALSE) standardGeneric("mapToGrid")
)

setGeneric("mapToRef", signature="grid",
    function(major, minor, grid, linear=FALSE) standardGeneric("mapToRef")
)

### Return an integer matrix with 1 column per dimension.
.normarg_Mindex <- function(Mindex, ndim, what="'Mindex'")
{
    if (!is.numeric(Mindex))
        stop(wmsg(what, " must be a numeric vector or matrix"))
    if (is.matrix(Mindex)) {
        if (ncol(Mindex) != ndim)
            stop(wmsg(what, " must have one column (or element ",
                      "if a vector) per dimension"))
    } else {
        if (is.array(Mindex))
            stop(wmsg(what, " must be a numeric vector or matrix"))
        if (length(Mindex) != ndim)
            stop(wmsg(what, " must have one element (or column ",
                      "if a matrix) per dimension"))
        Mindex <- matrix(Mindex, nrow=1)
    }
    if (storage.mode(Mindex) != "integer")
        storage.mode(Mindex) <- "integer"
    Mindex
}

### Return a numeric vector.
.normarg_ind <- function(ind, what="'ind'")
{
    if (!is.numeric(ind))
        stop(wmsg(what, " must be a numeric vector"))
    if (is.matrix(ind))
        stop(wmsg(what, " cannot be a matrix"))
    if (suppressWarnings(min(ind, na.rm=TRUE)) < 1)
        stop(wmsg(what, " must contain positive indices"))
    ind
}

.major_minor_as_list <- function(major, minor, grid, linear=FALSE)
{
    if (linear) {
        major <- Mindex2Lindex(major, dim(grid))
        ## 'as.integer=TRUE' is safe to use here because the viewports in
        ## an ArrayGrid object cannot be longer than 2^31 - 1 (this in turn
        ## translates in 'all(rowProds(dims(grid)) < .Machine$integer.max)').
        minor <- Mindex2Lindex(minor, dims(grid)[major, , drop=FALSE],
                               as.integer=TRUE)
    }
    list(major=major, minor=minor)
}

setMethod("mapToGrid", "ArbitraryArrayGrid",
    function(Mindex, grid, linear=FALSE)
    {
        if (!isTRUEorFALSE(linear))
            stop("'linear' must be TRUE or FALSE")
        ndim <- length(grid@tickmarks)
        Mindex <- .normarg_Mindex(Mindex, ndim)
        major <- lapply(seq_len(ndim),
            function(along) {
              1L + findInterval(Mindex[ , along], grid@tickmarks[[along]] + 1L)
            }
        )
        minor <- lapply(seq_len(ndim),
            function(along) {
                tm <- grid@tickmarks[[along]]
                tm_len <- length(tm)
                if (tm_len == 0L)
                    return(rep.int(NA_integer_, nrow(Mindex)))
                offset <- c(0L, tm[-tm_len])
                Mindex[ , along] - offset[major[[along]]]
            }
        )
        major <- do.call(cbind, major)
        minor <- do.call(cbind, minor)
        .major_minor_as_list(major, minor, grid, linear=linear)
    }
)

.normargs_major_minor <- function(major, minor, grid, linear)
{
    if (!isTRUEorFALSE(linear))
        stop("'linear' must be TRUE or FALSE")
    if (linear) {
        major <- .normarg_ind(major, what="when 'linear=TRUE', 'major'")
        minor <- .normarg_ind(minor, what="when 'linear=TRUE', 'minor'")
        if (length(major) != length(minor))
            stop(wmsg("when 'linear=TRUE', 'major' and 'minor' ",
                      "must have the same length"))
        minor <- Lindex2Mindex(minor, dims(grid)[major, , drop=FALSE])
        major <- Lindex2Mindex(major, dim(grid))
    } else {
        ndim <- length(refdim(grid))
        major <- .normarg_Mindex(major, ndim, what="'major'")
        minor <- .normarg_Mindex(minor, ndim, what="'minor'")
        if (nrow(major) != nrow(minor))
            stop(wmsg("'major' and 'minor' must have the same number of rows"))
    }
    list(major=major, minor=minor)
}

setMethod("mapToRef", "ArbitraryArrayGrid",
    function(major, minor, grid, linear=FALSE)
    {
        majmin <- .normargs_major_minor(major, minor, grid, linear)
        ans <- majmin$minor
        for (along in seq_len(ncol(ans))) {
            tm <- grid@tickmarks[[along]]
            tm_len <- length(tm)
            if (tm_len == 0L)
                next
            offset <- c(0L, tm[-tm_len])
            ans[ , along] <- ans[ , along] + offset[majmin$major[ , along]]
        }
        ans
    }
)

setMethod("mapToGrid", "RegularArrayGrid",
    function(Mindex, grid, linear=FALSE)
    {
        if (!isTRUEorFALSE(linear))
            stop("'linear' must be TRUE or FALSE")
        ndim <- length(grid@spacings)
        Mindex <- .normarg_Mindex(Mindex, ndim)
        d <- rep(grid@spacings, each=nrow(Mindex))
        Mindex0 <- Mindex - 1L  # 0-based indices
        major <- 1L + Mindex0 %/% d
        minor <- 1L + Mindex0 %% d
        .major_minor_as_list(major, minor, grid, linear=linear)
    }
)

setMethod("mapToRef", "RegularArrayGrid",
    function(major, minor, grid, linear=FALSE)
    {
        majmin <- .normargs_major_minor(major, minor, grid, linear)
        d <- rep(grid@spacings, each=nrow(majmin$major))
        (majmin$major - 1L) * d + majmin$minor
    }
)

