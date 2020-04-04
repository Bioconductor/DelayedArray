\name{ArrayGrid-class}
\docType{class}

\alias{class:ArrayViewport}
\alias{ArrayViewport-class}
\alias{ArrayViewport}

\alias{refdim}
\alias{refdim,ArrayViewport-method}
\alias{ranges,ArrayViewport-method}
\alias{start,ArrayViewport-method}
\alias{width,ArrayViewport-method}
\alias{end,ArrayViewport-method}
\alias{dim,ArrayViewport-method}

\alias{show,ArrayViewport-method}

\alias{makeNindexFromArrayViewport}

\alias{class:ArrayGrid}
\alias{ArrayGrid-class}
\alias{ArrayGrid}

\alias{class:ArbitraryArrayGrid}
\alias{ArbitraryArrayGrid-class}
\alias{ArbitraryArrayGrid}

\alias{class:RegularArrayGrid}
\alias{RegularArrayGrid-class}
\alias{RegularArrayGrid}

\alias{refdim,ArbitraryArrayGrid-method}
\alias{refdim,RegularArrayGrid-method}
\alias{dim,ArbitraryArrayGrid-method}
\alias{dim,RegularArrayGrid-method}
\alias{as.character.ArrayGrid}
\alias{as.character,ArrayGrid-method}
\alias{dims}
\alias{dims,ArrayGrid-method}
\alias{lengths,ArrayGrid-method}
\alias{maxlength}
\alias{maxlength,ANY-method}
\alias{maxlength,ArbitraryArrayGrid-method}
\alias{maxlength,RegularArrayGrid-method}
\alias{show,ArrayGrid-method}

\alias{aperm.ArbitraryArrayGrid}
\alias{aperm,ArbitraryArrayGrid-method}
\alias{aperm.RegularArrayGrid}
\alias{aperm,RegularArrayGrid-method}

\alias{downsample}
\alias{downsample,ArbitraryArrayGrid-method}
\alias{downsample,RegularArrayGrid-method}

\title{ArrayGrid and ArrayViewport objects}

\description{
  ArrayGrid and ArrayViewport objects are used internally to support
  block processing of array-like objects.
}

\usage{
## Constructor functions:
ArbitraryArrayGrid(tickmarks)
RegularArrayGrid(refdim, spacings=refdim)

downsample(x, ratio=1L)
}

\arguments{
  \item{tickmarks}{
    A list of integer vectors, one along each dimension of the reference
    array, representing the tickmarks along that dimension. Each integer
    vector must be sorted in ascending order. NAs or negative values are
    not allowed.
  }
  \item{refdim}{
    An integer vector containing the dimensions of the reference array.
  }
  \item{spacings}{
    An integer vector specifying the grid spacing along each dimension.
  }
  \item{x}{
    An ArrayGrid object.
  }
  \item{ratio}{
    An integer vector specifying the ratio of the downsampling along
    each dimension. Can be of length 1, in which case the same ratio is
    used along all the dimensions.
  }
}

\value{
  \itemize{
    \item For \code{ArbitraryArrayGrid()}: An ArbitraryArrayGrid instance.

    \item For \code{RegularArrayGrid()}: A RegularArrayGrid instance.

    \item For \code{downsample()}: An ArrayGrid object on the same reference
          array than \code{x}.
  }
}

\seealso{
  \itemize{
    \item \code{\link{read_block}}.

    \item \code{\link{blockApply}} for more information about block processing
          of array-like objects.

    \item \code{\link{mapToGrid}} for mapping reference array positions to
          grid positions and vice-versa.

    \item \code{\link{chunkGrid}}.

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. ArrayGrid OBJECTS
## ---------------------------------------------------------------------

## Create a regularly-spaced grid on top of a 3700 x 100 x 33 array:
grid1 <- RegularArrayGrid(c(3700, 100, 33), c(250, 100, 10))

## Dimensions of the reference array:
refdim(grid1)

## Number of grid elements along each dimension of the reference array:
dim(grid1)

## Total number of grid elements:
length(grid1)

## First element in the grid:
grid1[[1L]]             # same as grid1[[1L, 1L, 1L]]

## Last element in the grid:
grid1[[length(grid1)]]  # same as grid1[[15L, 1L, 4L]]

## Dimensions of the grid elements:
dims(grid1)             # one row per grid element

## Lengths of the grid elements:
lengths(grid1)          # same as rowProds(dims(grid1))
stopifnot(sum(lengths(grid1)) == prod(refdim(grid1)))

maxlength(grid1)        # does not need to compute lengths(grid1)) first
                        # so is more efficient than max(lengths(grid1))
stopifnot(maxlength(grid1) == max(lengths(grid1)))

## Create an arbitrary-spaced grid on top of a 15 x 9 matrix:
grid2 <- ArbitraryArrayGrid(list(c(2L, 7:10, 13L, 15L), c(5:6, 6L, 9L)))

refdim(grid2)
dim(grid2)
length(grid2)
grid2[[1L]]             # same as grid2[[1L, 1L]]
grid2[[length(grid2)]]  # same as grid2[[15L, 9L]]

dims(grid2)
lengths(grid2)
array(lengths(grid2), dim(grid2))  # display the grid element lengths in
                                   # an array of same shape as grid2

stopifnot(sum(lengths(grid2)) == prod(refdim(grid2)))

maxlength(grid2)        # does not need to compute lengths(grid2)) first
                        # so is more efficient than max(lengths(grid2))
stopifnot(maxlength(grid2) == max(lengths(grid2)))

## Max (i.e. highest) resolution grid:
Hgrid <- RegularArrayGrid(6:4, c(1, 1, 1))
Hgrid
dim(Hgrid)              # same as refdim(Hgrid)
stopifnot(identical(dim(Hgrid), refdim(Hgrid)))
stopifnot(all(lengths(Hgrid) == 1))

## Min (i.e. lowest) resolution grid:
Lgrid <- RegularArrayGrid(6:4, 6:4)
Lgrid
stopifnot(all(dim(Lgrid) == 1))
stopifnot(identical(dim(Lgrid[[1L]]), refdim(Lgrid)))
stopifnot(identical(dims(Lgrid), matrix(refdim(Lgrid), nrow=1)))

## ---------------------------------------------------------------------
## B. ArrayViewport OBJECTS
## ---------------------------------------------------------------------

## Grid elements are ArrayViewport objects:
grid1[[1L]]
class(grid1[[1L]])
grid1[[2L]]
grid1[[2L, 1L, 1L]]
grid1[[15L, 1L, 4L]]

## Construction of a standalong ArrayViewport object:
m0 <- matrix(1:30, ncol=5)
block_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m0), IRanges(c(3, 2), width=block_dim))
viewport1

dim(viewport1)     # 'block_dim'
length(viewport1)  # number of array elements in the viewport
ranges(viewport1)

## ---------------------------------------------------------------------
## C. GRIDS CAN BE TRANSPOSED
## ---------------------------------------------------------------------

tgrid2 <- t(grid2)
dim(tgrid2)
refdim(tgrid2)

## Use aperm() if the grid has more than 2 dimensions:
tgrid1 <- aperm(grid1)
dim(tgrid1)
refdim(tgrid1)

aperm(grid1, c(3, 1, 2))
aperm(grid1, c(1, 3, 2))
aperm(grid1, c(3, 1))     # some dimensions can be dropped
aperm(grid1, c(3, 2, 3))  # and some can be repeated

## ---------------------------------------------------------------------
## D. DOWNSAMPLING AN ArrayGrid OBJECT
## ---------------------------------------------------------------------
## The elements (ArrayViewport) of an ArrayGrid object can be replaced
## with bigger elements obtained by merging adjacent elements. How many
## adjacent elements to merge along each dimension is specified via the
## 'ratio' vector (one integer per dimension). We call this operation
## "downsampling. It can be seen as reducing the "resolution" of a grid
## by the specified ratio (if we think of the grid elements as pixels).
downsample(grid2, 2)
downsample(grid2, 3)
downsample(grid2, 4)

## Downsampling preserves the dimensions of the reference array:
stopifnot(identical(refdim(downsample(grid2, 2)), refdim(grid2)))
stopifnot(identical(refdim(downsample(grid2, 3)), refdim(grid2)))
stopifnot(identical(refdim(downsample(grid2, 4)), refdim(grid2)))

## A big enough ratio will eventually produce the coarsest possible grid
## i.e. a grid with a single grid element covering the entire reference
## array:
grid3 <- downsample(grid2, 7)
length(grid3)
grid3[[1L]]
stopifnot(identical(dim(grid3[[1L]]), refdim(grid3)))

## Downsampling by a ratio of 1 is a no-op:
stopifnot(identical(downsample(grid2, 1), grid2))

## Using one ratio per dimension:
downsample(grid2, c(2, 1))

## Downsample a max resolution grid:
refdim <- c(45, 16, 20)
grid4 <- RegularArrayGrid(refdim, c(1, 1, 1))
ratio <- c(6, 1, 3)
stopifnot(identical(
    downsample(grid4, ratio),
    RegularArrayGrid(refdim, ratio)
))
}
\keyword{classes}
\keyword{methods}
