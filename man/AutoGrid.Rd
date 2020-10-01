\name{AutoGrid}

\alias{AutoGrid}

\alias{defaultAutoGrid}
\alias{blockGrid}
\alias{rowAutoGrid}
\alias{rowGrid}
\alias{colAutoGrid}
\alias{colGrid}
\alias{getAutoGridMaker}
\alias{setAutoGridMaker}

\title{Generate automatic grids to use for block processing
       of array-like objects}

\description{
  \code{defaultAutoGrid()} is the default \emph{automatic grid maker}.
  It creates a grid that is suitable for block processing of the array-like
  object passed to it.

  \code{rowAutoGrid()} and \code{colAutoGrid()} are additional
  \emph{automatic grid makers}, specific to the 2-dimensional case.
  They can be used to define blocks of full rows or full columns.
}

\usage{
defaultAutoGrid(x, block.length=NULL, chunk.grid=NULL, block.shape=NULL)

rowAutoGrid(x, nrow=NULL, block.length=NULL)
colAutoGrid(x, ncol=NULL, block.length=NULL)

## Replace default automatic grid maker with user-defined one:
getAutoGridMaker()
setAutoGridMaker(GRIDMAKER="defaultAutoGrid")
}

\arguments{
  \item{x}{
    An array-like or matrix-like object for \code{defaultAutoGrid}.

    A matrix-like object for \code{rowAutoGrid} and \code{colAutoGrid}.
  }
  \item{block.length}{
    The length of the blocks i.e. the number of array elements per block.
    By default the automatic block length (returned by
    \code{getAutoBlockLength(type(x))}) is used.
    Depending on how much memory is available on your machine, you might
    want to increase (or decrease) the automatic block length by adjusting
    the automatic block size with \code{setAutoBlockSize()}.
  }
  \item{chunk.grid}{
    The grid of physical chunks.
    By default \code{\link{chunkGrid}(x)} is used.
  }
  \item{block.shape}{
    A string specifying the shape of the blocks.
    See \code{\link{makeCappedVolumeBox}} for a description of the
    supported shapes.
    By default \code{getAutoBlockShape()} is used.
  }
  \item{nrow}{
    The number of rows of the blocks. The bottommost blocks might have less.
    See examples below.
  }
  \item{ncol}{
    The number of columns of the blocks. The rightmost blocks might have less.
    See examples below.
  }
  \item{GRIDMAKER}{
    The function to use as \emph{automatic grid maker}, that is, the
    function that will be used by \code{\link{blockApply}()} and
    \code{\link{blockReduce}()} to make a grid when no grid is supplied
    via their \code{grid} argument.
    The function will be called on array-like object \code{x} and must
    return an \link{ArrayGrid} object, say \code{grid}, that is compatible
    with \code{x} i.e. such that \code{refdim(grid)} is identical to
    \code{dim(x)}.

    \code{GRIDMAKER} can be specified as a function or as a single string
    naming a function. It can be a user-defined function or a pre-defined
    grid maker like \code{defaultAutoGrid}, \code{rowAutoGrid}, or
    \code{colAutoGrid}.

    The \emph{automatic grid maker} is set to \code{defaultAutoGrid} at
    package startup and can be reset anytime to this value by calling
    \code{setAutoGridMaker()} with no argument.
  }
}

\details{
  By default, primary block processing functions \code{\link{blockApply}()}
  and \code{\link{blockReduce}()} use the grid returned by
  \code{defaultAutoGrid(x)} to process array-like object \code{x}
  block by block. This can be changed with \code{setAutoGridMaker()}.
}

\value{
  \code{defaultAutoGrid}: An \link{ArrayGrid} object on reference array
  \code{x}. The grid elements define the blocks that will be used to
  process \code{x} by block. The grid is \emph{optimal} in the sense that:
  \enumerate{
    \item It's \emph{compatible} with the grid of physical chunks a.k.a.
          \emph{chunk grid}. This means that, when the chunk grid is known
          (i.e. when \code{\link{chunkGrid}(x)} is not NULL or
          \code{chunk.grid} is supplied), every block in the grid contains
          one or more \emph{full} chunks. In other words, chunks never cross
          block boundaries.
    \item Its \emph{resolution} is such that the blocks have a length
          that is as close as possibe to (but does not exceed)
          \code{block.length}. An exception is made when some chunks
          already have a length that is >= \code{block.length}, in which
          case the returned grid is the same as the chunk grid.
  }
  Note that the returned grid is regular (i.e. is a \link{RegularArrayGrid}
  object) unless the chunk grid is not regular (i.e. is an
  \link{ArbitraryArrayGrid} object).

  \code{rowAutoGrid}: A \link{RegularArrayGrid} object on reference array
  \code{x} where the grid elements define blocks made of full rows of \code{x}.

  \code{colAutoGrid}: A \link{RegularArrayGrid} object on reference array
  \code{x} where the grid elements define blocks made of full columns
  of \code{x}.
}

\seealso{
  \itemize{
    \item \code{\link{setAutoBlockSize}} and \code{\link{setAutoBlockShape}}
          to control the geometry of automatic blocks.

    \item \code{\link{blockApply}} for convenient block processing of an
          array-like object.

    \item \link{ArrayGrid} objects.

    \item The \code{\link{makeCappedVolumeBox}} utility to make
          \emph{capped volume boxes}.

    \item \code{\link{chunkGrid}}.

    \item \code{\link{read_block}} and \code{\link{write_block}}.
  }
}

\examples{
## ---------------------------------------------------------------------
## A VERSION OF sum() THAT USES BLOCK PROCESSING
## ---------------------------------------------------------------------

block_sum <- function(a, grid) {
  sums <- lapply(grid, function(viewport) sum(read_block(a, viewport)))
  sum(unlist(sums))
}

## On an ordinary matrix:
m <- matrix(runif(600), ncol=12)
m_grid <- defaultAutoGrid(m, block.length=120)
sum1 <- block_sum(m, m_grid)
sum1

## On a DelayedArray object:
library(HDF5Array)
M <- as(m, "HDF5Array")
sum2 <- block_sum(M, m_grid)
sum2

sum3 <- block_sum(M, colAutoGrid(M, block.length=120))
sum3

sum4 <- block_sum(M, rowAutoGrid(M, block.length=80))
sum4

## Sanity checks:
sum0 <- sum(m)
stopifnot(identical(sum1, sum0))
stopifnot(identical(sum2, sum0))
stopifnot(identical(sum3, sum0))
stopifnot(identical(sum4, sum0))

## ---------------------------------------------------------------------
## defaultAutoGrid()
## ---------------------------------------------------------------------
grid <- defaultAutoGrid(m, block.length=120)
grid
as.list(grid)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- defaultAutoGrid(m, block.length=120,
                           block.shape="first-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- defaultAutoGrid(m, block.length=120,
                           block.shape="last-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

defaultAutoGrid(m, block.length=100)
defaultAutoGrid(m, block.length=75)
defaultAutoGrid(m, block.length=25)
defaultAutoGrid(m, block.length=20)
defaultAutoGrid(m, block.length=10)

## ---------------------------------------------------------------------
## rowAutoGrid() AND colAutoGrid()
## ---------------------------------------------------------------------
rowAutoGrid(m, nrow=10)  # 5 blocks of 10 rows each
rowAutoGrid(m, nrow=15)  # 3 blocks of 15 rows each plus 1 block of 5 rows
colAutoGrid(m, ncol=5)   # 2 blocks of 5 cols each plus 1 block of 2 cols

## See '?write_block' for an advanced example of user-implemented
## block processing using colAutoGrid() and a realization sink.

## ---------------------------------------------------------------------
## REPLACE DEFAULT AUTOMATIC GRID MAKER WITH USER-DEFINED ONE
## ---------------------------------------------------------------------
getAutoGridMaker()
setAutoGridMaker(function(x) colAutoGrid(x, ncol=5))
getAutoGridMaker()

blockApply(m, function(block) currentViewport())

## Reset automatic grid maker to factory settings:
setAutoGridMaker()
}
\keyword{utilities}
