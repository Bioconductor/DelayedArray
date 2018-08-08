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
## A. blockGrid()
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
## B. rowGrid() AND colGrid()
## ---------------------------------------------------------------------
rowGrid(m, nrow=10)  # 5 blocks of 10 rows each
rowGrid(m, nrow=15)  # 3 blocks of 15 rows each plus 1 block of 5 rows
colGrid(m, ncol=5)   # 2 blocks of 5 cols each plus 1 block of 2 cols

## ---------------------------------------------------------------------
## C. AN ADVANCED EXAMPLE OF USER-CONTROLLED BLOCK PROCESSING USING
##    colGrid() AND A REALIZATION SINK
## ---------------------------------------------------------------------
## Say we have 2 matrices with the same number of columns. Each column
## represents a biological sample:
library(HDF5Array)
R <- as(matrix(runif(75000), ncol=1000), "HDF5Array")   # 75 rows
G <- as(matrix(runif(250000), ncol=1000), "HDF5Array")  # 250 rows

## Say we want to compute the matrix U obtained by applying the same
## binary functions FUN() to all samples i.e. U is defined as:
##
##   U[ , j] <- FUN(R[ , j], G[ , j]) for 1 <= j <= 1000
##
## Note that FUN() should return a vector of constant length, say 200,
## so U will be a 200x1000 matrix. A naive implementation would be:
##
##   pFUN <- function(r, g) {
##       stopifnot(ncol(r) == ncol(g))  # sanity check
##       sapply(seq_len(ncol(r)), function(j) FUN(r[ , j], g[ , j]))
##   }
##
## But because U is going to be too big to fit in memory, we can't
## just do pFUN(R, G). So we want to compute U block by block and
## write the blocks to disk as we go. The blocks will be made of full
## columns. Also since we need to walk on 2 matrices at the same time
## (R and G), we can't use blockApply() or blockReduce() so we'll use
## a "for" loop.

## Before we can write the "for", we need 4 things:

## 1) Two grids of blocks, one on R and one on G. The blocks in the
##    2 grids must contain the same number of columns. We arbitrarily
##    choose to use blocks of 150 columns:
R_grid <- colGrid(R, ncol=150)
G_grid <- colGrid(G, ncol=150)

## 2) The function pFUN(). It will take 2 blocks as input, 1 from R
##    and 1 from G, apply FUN() to all the samples in the blocks,
##    and return a matrix with one columns per sample:
pFUN <- function(r, g) {
    stopifnot(ncol(r) == ncol(g))  # sanity check
    ## Return a matrix with 200 rows with random values. Completely
    ## artificial sorry. A realistic example would actually need to
    ## apply the same binary function to r[ ,j] and g[ , j] for
    ## 1 <= j <= ncol(r).
    matrix(runif(200 * ncol(r)), nrow=200)
}

## 3) A RealizationSink object where to write the matrices returned
##    by pFUN() as we go. Note that instead of creating a realization
##    sink by calling a backend-specific sink constructor (e.g.
##    HDF5Array:::HDF5RealizationSink), we use the backend-agnostic
##    constructor DelayedArray:::RealizationSink() and set the current
##    realization backend to HDF5:

setRealizationBackend("HDF5Array")
U_sink <- DelayedArray:::RealizationSink(c(200L, 1000L))

## 4) Finally, we create a grid on U_sink with blocks that contain the
##    same number of columns as the corresponding blocks in R and G:

U_grid <- colGrid(U_sink, ncol=150)

## Note that the 3 grids should have the same number of blocks:
stopifnot(length(U_grid) == length(R_grid))
stopifnot(length(U_grid) == length(G_grid))

## Now we can procede. We write a loop where we walk on R and G at
## the same time, block by block, apply pFUN(), and write the output
## of pFUN() to U_sink:
for (b in seq_along(U_grid)) {
    R_block <- read_block(R, R_grid[[b]])
    G_block <- read_block(G, G_grid[[b]])
    U_block <- pFUN(R_block, G_block)
    write_block(U_sink, U_grid[[b]], U_block)
}

close(U_sink)
U <- as(U_sink, "DelayedArray")

## A note about parallelization: even though concurrent block reading
## from the same object is supported, concurrent writing to a sink is
## not supported yet. So the above code cannot be parallelized at the
## moment.
}
\keyword{utilities}
