### =========================================================================
### Map reference array positions to grid positions and vice-versa
### -------------------------------------------------------------------------


### 'aind' must be a numeric vector or matrix (a vector is treated
### like a 1-row matrix).
setGeneric("mapToGrid", signature="grid",
    function(aind, grid, linear=FALSE) standardGeneric("mapToGrid")
)

setGeneric("mapToRef", signature="grid",
    function(major, minor, grid, linear=FALSE) standardGeneric("mapToRef")
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
                1L + findInterval(aind[ , along], grid@tickmarks[[along]] + 1L)
            }
        )
        minor <- lapply(seq_len(ndim),
            function(along) {
                tm <- grid@tickmarks[[along]]
                tm_len <- length(tm)
                if (tm_len == 0L)
                    return(rep.int(NA_integer_, nrow(aind)))
                offset <- c(0L, tm[-tm_len])
                aind[ , along] - offset[major[[along]]]
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
        major <- normarg_ind(major, what="when 'linear=TRUE', 'major'")
        minor <- normarg_ind(minor, what="when 'linear=TRUE', 'minor'")
        if (length(major) != length(minor))
            stop(wmsg("when 'linear=TRUE', 'major' and 'minor' ",
                      "must have the same length"))
        minor <- arrayInd2(minor, dims(grid)[major, , drop=FALSE])
        major <- arrayInd2(major, dim(grid))
    } else {
        ndim <- length(refdim(grid))
        major <- normarg_aind(major, ndim, what="'major'")
        minor <- normarg_aind(minor, ndim, what="'minor'")
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

setMethod("mapToRef", "RegularArrayGrid",
    function(major, minor, grid, linear=FALSE)
    {
        majmin <- .normargs_major_minor(major, minor, grid, linear)
        d <- rep(grid@spacings, each=nrow(majmin$major))
        (majmin$major - 1L) * d + majmin$minor
    }
)

