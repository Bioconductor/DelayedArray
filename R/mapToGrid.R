### =========================================================================
### mapToGrid()
### -------------------------------------------------------------------------


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

