\name{blockGrid}

\alias{blockGrid}
\alias{rowGrid}
\alias{colGrid}

\title{Define a grid suitable for block processing of an array-like object}

\description{
  \code{blockGrid()} is the primary utility function to use to define a grid
  that is suitable for block processing of an array-like object.

  \code{rowGrid()} and \code{colGrid()} are additional functions, specific
  to the 2-dimensional case, that can be used to define blocks of full rows
  or full columns.
}

\usage{
blockGrid(x, block.maxlength=NULL, chunk.grid=NULL,
             block.shape=c("hypercube",
                           "scale",
                           "first-dim-grows-first",
                           "last-dim-grows-first"))

rowGrid(x, nrow=NULL, block.maxlength=NULL)
colGrid(x, ncol=NULL, block.maxlength=NULL)
}

\arguments{
  \item{x}{
    An array-like or matrix-like object for \code{blockGrid}.

    A matrix-like object for \code{rowGrid} and \code{colGrid}.
  }
  \item{block.maxlength}{
    The maximum length of the blocks.
    By default \code{getOption("DelayedArray.block.size") / sizeof(type(x))}
    is used. In a fresh session (i.e. when the \pkg{DelayedArray} package is
    loaded), global option \code{DelayedArray.block.size} is set to 45e6
    or 45 Mb (this is internal constant
    \code{DelayedArray:::DEFAULT_BLOCK_SIZE}).
    Depending on how much memory is available on your machine, you can
    increase or decrease the default block size with something like
    \code{options(DelayedArray.block.size=200e6)}.
  }
  \item{chunk.grid}{
    The grid of physical chunks.
    By default \code{\link{chunkGrid}(x)} is used.
  }
  \item{block.shape}{
    A string specifying the shape of the blocks.
    See \code{\link{makeCappedVolumeBox}} for a description of the supported
    shapes.
  }
  \item{nrow}{
    The number of rows of the blocks. The bottommost blocks might have less.
    See examples below.
  }
  \item{ncol}{
    The number of columns of the blocks. The rightmost blocks might have less.
    See examples below.
  }
}

\details{
  By default, primary block processing functions \code{\link{blockApply}()}
  and \code{\link{blockReduce}()} use the grid returned by \code{blockGrid(x)}
  to process array-like object \code{x} block by block.
}

\value{
  \code{blockGrid}: An \link{ArrayGrid} object on reference array \code{x}.
  The grid elements define the blocks that will be used to process \code{x}
  by block. The grid is \emph{optimal} in the sense that:
  \enumerate{
    \item It's \emph{compatible} with the grid of physical chunks a.k.a.
          \emph{chunk grid}. This means that, when the chunk grid is known
          (i.e. when \code{\link{chunkGrid}(x)} is not NULL or
          \code{chunk.grid} is supplied), every block in the grid contains
          one or more \emph{full} chunks. In other words, chunks never cross
          block boundaries.
    \item Its \emph{resolution} is such that the blocks have a length
          that is as close as possibe to (but does not exceed)
          \code{block.maxlength}. An exception is made when some chunks
          already have a length that is >= \code{block.maxlength}, in which
          case the returned grid is the same as the chunk grid.
  }
  Note that the returned grid is regular (i.e. is a \link{RegularArrayGrid}
  object) unless the chunk grid is not regular (i.e. is an
  \link{ArbitraryArrayGrid} object).

  \code{rowGrid}: A \link{RegularArrayGrid} object on reference array \code{x}
  where the grid elements define blocks made of full rows of \code{x}.

  \code{colGrid}: A \link{RegularArrayGrid} object on reference array \code{x}
  where the grid elements define blocks made of full columns of \code{x}.
}

\seealso{
  \itemize{
    \item \code{\link{blockApply}} and family to process an array-like
          object block by block.

    \item \link{ArrayGrid} objects.

    \item The \code{\link{makeCappedVolumeBox}} utility to make
          \emph{capped volume boxes}.

    \item \code{\link{chunkGrid}}.
  }
}

\examples{
m <- matrix(runif(600), ncol=12)

## ---------------------------------------------------------------------
## blockGrid()
## ---------------------------------------------------------------------
grid <- blockGrid(m, block.maxlength=120)
grid
as.list(grid)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- blockGrid(m, block.maxlength=120,
                     block.shape="first-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- blockGrid(m, block.maxlength=120,
                     block.shape="last-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

blockGrid(m, block.maxlength=100)
blockGrid(m, block.maxlength=75)
blockGrid(m, block.maxlength=25)
blockGrid(m, block.maxlength=20)
blockGrid(m, block.maxlength=10)

## ---------------------------------------------------------------------
## rowGrid() AND colGrid()
## ---------------------------------------------------------------------
rowGrid(m, nrow=10)  # 5 blocks of 10 rows each
rowGrid(m, nrow=15)  # 3 blocks of 15 rows each plus 1 block of 5 rows
colGrid(m, ncol=5)   # 2 blocks of 5 cols each plus 1 block of 2 cols
}
\keyword{utilities}
