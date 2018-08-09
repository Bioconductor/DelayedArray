\name{blockGrid}

\alias{blockGrid}
\alias{rowGrid}
\alias{colGrid}

\alias{getDefaultBlockSize}
\alias{setDefaultBlockSize}
\alias{getDefaultBlockLength}
\alias{getDefaultBlockShape}
\alias{setDefaultBlockShape}

\title{Define grids to use in the context of block processing}

\description{
  \code{blockGrid()} is the primary utility function to use to define a grid
  that is suitable for block processing of an array-like object.

  \code{rowGrid()} and \code{colGrid()} are additional functions, specific
  to the 2-dimensional case. They can be used to define blocks of full rows
  or full columns.

  A family of utilities is provided to control the default block size (or
  length) and shape.
}

\usage{
## Define grids to use in the context of block processing:

blockGrid(x, block.length=NULL, chunk.grid=NULL, block.shape=NULL)

rowGrid(x, nrow=NULL, block.length=NULL)
colGrid(x, ncol=NULL, block.length=NULL)

## Control the default block size (or length) and shape:

getDefaultBlockSize()
setDefaultBlockSize(size=1e8)

getDefaultBlockLength(type)

getDefaultBlockShape()
setDefaultBlockShape(shape=c("hypercube",
                             "scale",
                             "first-dim-grows-first",
                             "last-dim-grows-first"))
}

\arguments{
  \item{x}{
    An array-like or matrix-like object for \code{blockGrid}.

    A matrix-like object for \code{rowGrid} and \code{colGrid}.
  }
  \item{block.length}{
    The length of the blocks i.e. the number of array elements per block.
    By default \code{getDefaultBlockLength(type(x))} is used.
    Depending on how much memory is available on your machine, you might
    want to increase (or decrease) the default block size with something
    like \code{setDefaultBlockSize(5e8)}.
  }
  \item{chunk.grid}{
    The grid of physical chunks.
    By default \code{\link{chunkGrid}(x)} is used.
  }
  \item{block.shape}{
    A string specifying the shape of the blocks.
    See \code{\link{makeCappedVolumeBox}} for a description of the
    supported shapes.
    By default \code{getDefaultBlockShape()} is used.
  }
  \item{nrow}{
    The number of rows of the blocks. The bottommost blocks might have less.
    See examples below.
  }
  \item{ncol}{
    The number of columns of the blocks. The rightmost blocks might have less.
    See examples below.
  }
  \item{size}{
    The size of the blocks in bytes. Note that, except when the type of the
    array data is \code{"character"}, the size of a block is its length
    multiplied by the size of an array element. For example, a block of
    500x1000x500 doubles has a length of 250 million elements and a size
    of 2 Gb (each double occupies 8 bytes of memory).

    The default block size is set to 100 Mb at package startup and can
    be reset anytime to this value by calling \code{setDefaultBlockSize()}
    with no argument.
  }
  \item{type}{
    A string specifying the type of the array data.
  }
  \item{shape}{
    A string specifying the shape of the blocks.
    See \code{\link{makeCappedVolumeBox}} for a description of the
    supported shapes.

    The default block shape is set to \code{"hypercube"} at package startup
    and can be reset anytime to this value by calling
    \code{setDefaultBlockShape()} with no argument.
  }
}

\details{
  By default, primary block processing functions \code{\link{blockApply}()}
  and \code{\link{blockReduce}()} use the grid returned by \code{blockGrid(x)}
  to process array-like object \code{x} block by block. This can be changed
  with \code{\link{setDefaultGridMaker}()}.
  See \code{?\link{setDefaultGridMaker}} for more information.
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
          \code{block.length}. An exception is made when some chunks
          already have a length that is >= \code{block.length}, in which
          case the returned grid is the same as the chunk grid.
  }
  Note that the returned grid is regular (i.e. is a \link{RegularArrayGrid}
  object) unless the chunk grid is not regular (i.e. is an
  \link{ArbitraryArrayGrid} object).

  \code{rowGrid}: A \link{RegularArrayGrid} object on reference array \code{x}
  where the grid elements define blocks made of full rows of \code{x}.

  \code{colGrid}: A \link{RegularArrayGrid} object on reference array \code{x}
  where the grid elements define blocks made of full columns of \code{x}.

  \code{getDefaultBlockSize}: The current default block size in bytes
  as a single numeric value.

  \code{setDefaultBlockSize}: The new default block size in bytes as an
  invisible single numeric value.

  \code{getDefaultBlockLength}: The default block length as a single integer
  value.

  \code{getDefaultBlockShape}: The current default block shape as a
  single string.

  \code{setDefaultBlockShape}: The new default block shape as an invisible
  single string.
}

\seealso{
  \itemize{
    \item \code{\link{blockApply}} and family to process an array-like
          object block by block.

    \item \link{ArrayGrid} objects.

    \item The \code{\link{makeCappedVolumeBox}} utility to make
          \emph{capped volume boxes}.

    \item \code{\link{chunkGrid}}.

    \item Advanced users: \link{RealizationSink} objects for writing an 
          array-like object block by block to disk (or to memory).
  }
}

\examples{
m <- matrix(runif(600), ncol=12)

## ---------------------------------------------------------------------
## A. blockGrid()
## ---------------------------------------------------------------------
grid <- blockGrid(m, block.length=120)
grid
as.list(grid)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- blockGrid(m, block.length=120,
                     block.shape="first-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

grid <- blockGrid(m, block.length=120,
                     block.shape="last-dim-grows-first")
grid
table(lengths(grid))
stopifnot(maxlength(grid) <= 120)

blockGrid(m, block.length=100)
blockGrid(m, block.length=75)
blockGrid(m, block.length=25)
blockGrid(m, block.length=20)
blockGrid(m, block.length=10)

## ---------------------------------------------------------------------
## B. rowGrid() AND colGrid()
## ---------------------------------------------------------------------
rowGrid(m, nrow=10)  # 5 blocks of 10 rows each
rowGrid(m, nrow=15)  # 3 blocks of 15 rows each plus 1 block of 5 rows
colGrid(m, ncol=5)   # 2 blocks of 5 cols each plus 1 block of 2 cols

## See ?RealizationSink for an advanced example of user-controlled block
## processing using colGrid() and a realization sink.

## ---------------------------------------------------------------------
## C. CONTROL THE DEFAULT BLOCK SIZE (OR LENGTH) AND SHAPE
## ---------------------------------------------------------------------
getDefaultBlockSize()

getDefaultBlockLength("double")
getDefaultBlockLength("integer")
getDefaultBlockLength("logical")
getDefaultBlockLength("raw")

setDefaultBlockSize(140)
getDefaultBlockLength(type(m))
blockGrid(m)
lengths(blockGrid(m))
dims(blockGrid(m))

getDefaultBlockShape()
setDefaultBlockShape("scale")
blockGrid(m)
lengths(blockGrid(m))
dims(blockGrid(m))

## Reset default block size and shape to factory settings:
setDefaultBlockSize()
setDefaultBlockShape()
}
\keyword{utilities}
