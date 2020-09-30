\name{blockApply}

\alias{block processing}
\alias{block_processing}

\alias{blockApply}
\alias{blockReduce}

\alias{effectiveGrid}
\alias{currentBlockId}
\alias{currentViewport}

\alias{getAutoBPPARAM}
\alias{setAutoBPPARAM}

\title{Convenient block processing of an array-like object}

\description{
  Conveniently walk on the blocks of an array-like object and process them.
}

\usage{
blockApply(x, FUN, ..., grid=NULL, as.sparse=FALSE,
              BPPARAM=getAutoBPPARAM())

blockReduce(FUN, x, init, BREAKIF=NULL, grid=NULL, as.sparse=FALSE)

## Retrieve grid context for the current block:
effectiveGrid(envir=parent.frame(2))
currentBlockId(envir=parent.frame(2))
currentViewport(envir=parent.frame(2))

## Get/set automatic parallel back-end:
getAutoBPPARAM()
setAutoBPPARAM(BPPARAM=NULL)
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{FUN}{
    The function to be applied to each block of \code{x}. It must be able
    to accept as input any of the blocks of \code{x}.

    For \code{blockReduce()}, \code{init <- FUN(block, init)} will be
    performed on each block so \code{FUN} must take at least two arguments
    (typically \code{block} and \code{init} but the names can differ)
    and should normally return a value of the same type as its 2nd
    argument (\code{init}).

    IMPORTANT: If \code{as.sparse} is set to \code{FALSE}, all blocks will
    be passed to \code{FUN} as ordinary arrays. If it's set to \code{TRUE},
    they will be passed as \link{SparseArraySeed} objects.
    If it's set to \code{NA}, then \code{is_sparse(x)} determines how they
    will be passed.
  }
  \item{...}{
    Optional arguments to \code{FUN}.
  }
  \item{grid}{
    An \link{ArrayGrid} object compatible with the geometry of \code{x}.
    It defines the blocks to walk on. If not specified, an automatic grid
    is used. By default \code{\link{defaultAutoGrid}(x)} is called to
    generate the automatic grid. The \emph{automatic grid maker} can be
    changed with \code{\link{setAutoGridMaker}()}.
    See \code{?\link{setAutoGridMaker}} for more information.
  }
  \item{as.sparse}{
    Passed to the internal calls to \code{read_block}.
    See \code{?\link{read_block}} for more information.
  }
  \item{BPPARAM}{
    A \code{NULL}, in which case blocks are processed sequentially, or
    a \link[BiocParallel]{BiocParallelParam} instance (from the
    \pkg{BiocParallel} package), in which case they are processed in
    parallel. The specific \link[BiocParallel]{BiocParallelParam}
    instance determines the parallel back-end to use.
    See \code{?\link[BiocParallel]{BiocParallelParam}} in the
    \pkg{BiocParallel} package for more information about back-ends.
  }
  \item{init}{
    The value to passe to the first call to \code{FUN(block, init)} when
    \code{Reduce()} starts the walk. Note that \code{Reduce()} always
    operates sequentially.
  }
  \item{BREAKIF}{
    An optional callback function that detects a break condition.
    Must return \code{TRUE} or \code{FALSE}.
    At each iteration \code{blockReduce()} will call it on the
    result of \code{init <- FUN(block, init)} and exit the walk
    if \code{BREAKIF(init)} returned \code{TRUE}.
  }
  \item{envir}{
    Do not use (unless you know what you are doing).
  }
}

\details{
  \code{effectiveGrid()}, \code{currentBlockId()}, and \code{currentViewport()}
  can only be called (with no arguments!) from within the callback functions
  passed to \code{blockApply()} and \code{blockReduce()}, that is, from within
  functions \code{FUN} and \code{BREAKIF}.
  They return the grid context for the block being currently processed.
  See Value section below for the details.
}

\value{
  For \code{blockApply()}, a list with one list element per block visited.

  For \code{Reduce()}, the result of the last call to \code{FUN(block, init)}.

  For \code{effectiveGrid()}, the grid (\link{ArrayGrid} object) being
  effectively used.

  For \code{currentBlockId()}, the id (a.k.a. rank) of the current block.

  For \code{currentViewport()}, the viewport (\link{ArrayViewport} object)
  of the current block.
}

\seealso{
  \itemize{
    \item \code{\link{defaultAutoGrid}} and family to define grids to use
          in the context of block processing of array-like objects.

    \item \link{ArrayGrid} and \link{ArrayViewport} objects.

    \item \code{\link{read_block}}.

    \item \code{\link[BiocParallel]{MulticoreParam}},
          \code{\link[BiocParallel]{SnowParam}}, and
          \code{\link[BiocParallel]{bpparam}}, from the \pkg{BiocParallel}
          package.

    \item \link{DelayedArray} objects.

    \item Advanced users: \link{RealizationSink} objects for writing
          an array-like object block by block to disk (or to memory).
  }
}

\examples{
m <- matrix(1:60, nrow=10)
m_grid <- defaultAutoGrid(m, block.length=16, block.shape="hypercube")

## ---------------------------------------------------------------------
## blockApply()
## ---------------------------------------------------------------------
blockApply(m, identity, grid=m_grid)
blockApply(m, sum, grid=m_grid)

blockApply(m, function(block) {block + currentBlockId()*1e3}, grid=m_grid)
blockApply(m, function(block) currentViewport(), grid=m_grid)
blockApply(m, dim, grid=m_grid)

## The grid does not need to be regularly spaced:
a <- array(runif(8000), dim=c(25, 40, 8))
a_tickmarks <- list(c(7L, 15L, 25L), c(14L, 22L, 40L), c(2L, 8L))
a_grid <- ArbitraryArrayGrid(a_tickmarks)
a_grid
blockApply(a, function(block) sum(log(block + 0.5)), grid=a_grid)

## See block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
blockApply(m, function(block) sum(log(block + 0.5)), grid=m_grid)

## Use parallel evaluation:
library(BiocParallel)
if (.Platform$OS.type != "windows") {
    BPPARAM <- MulticoreParam(workers=4)
} else {
    ## MulticoreParam() is not supported on Windows so we use
    ## SnowParam() on this platform.
    BPPARAM <- SnowParam(4)
}
blockApply(m, function(block) sum(log(block + 0.5)), grid=m_grid,
              BPPARAM=BPPARAM)
## Note that blocks can be visited in any order!

## ---------------------------------------------------------------------
## blockReduce()
## ---------------------------------------------------------------------
FUN <- function(block, init) anyNA(block) || init
blockReduce(FUN, m, init=FALSE, grid=m_grid)

m[10, 1] <- NA
blockReduce(FUN, m, init=FALSE, grid=m_grid)

## With early bailout:
blockReduce(FUN, m, init=FALSE, BREAKIF=identity, grid=m_grid)

## Note that this is how the anyNA() method for DelayedArray object is
## implemented.

DelayedArray:::set_verbose_block_processing(FALSE)
}
\keyword{methods}
