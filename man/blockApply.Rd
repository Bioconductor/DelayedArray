\name{blockApply}

\alias{block processing}
\alias{block_processing}

\alias{setAutoBPPARAM}
\alias{getAutoBPPARAM}

\alias{set_grid_context}
\alias{effectiveGrid}
\alias{currentBlockId}
\alias{currentViewport}

\alias{viewportApply}
\alias{blockApply}
\alias{viewportReduce}
\alias{blockReduce}

\title{Convenient block processing of an array-like object}

\description{
  Conveniently walk on the blocks of an array-like object and process them.
}

\usage{
## Main looping functions:

blockApply(x, FUN, ..., grid=NULL, as.sparse=FALSE,
           BPPARAM=getAutoBPPARAM(), verbose=NA)

blockReduce(FUN, x, init, ..., BREAKIF=NULL, grid=NULL, as.sparse=FALSE,
            verbose=NA)

## Lower-level looping functions:
viewportApply(grid, FUN, ..., BPPARAM=getAutoBPPARAM(), verbose=NA)
viewportReduce(FUN, grid, init, ..., BREAKIF=NULL, verbose=NA)

## Retrieve grid context for the current block/viewport:
effectiveGrid(envir=parent.frame(2))
currentBlockId(envir=parent.frame(2))
currentViewport(envir=parent.frame(2))

## Get/set automatic parallel back-end:
getAutoBPPARAM()
setAutoBPPARAM(BPPARAM=NULL)

## For testing/debugging callback functions:
set_grid_context(effective_grid, current_block_id, envir=parent.frame(1))
}

\arguments{
  \item{x}{
    An array-like object, typically a \link{DelayedArray} object
    or derivative.
  }
  \item{FUN}{
    For \code{blockApply} and \code{blockReduce}, \code{FUN} is the callback
    function to apply to each block of \code{x}. It must be able to accept as
    input any of the blocks of \code{x}.

    IMPORTANT: If \code{as.sparse} is set to \code{FALSE}, all blocks will
    be passed to \code{FUN} as ordinary arrays. If it's set to \code{TRUE},
    they will be passed as \link{SparseArraySeed} objects.
    If it's set to \code{NA}, then \code{is_sparse(x)} determines how they
    will be passed to \code{FUN}.

    For \code{viewportApply()} and \code{viewportReduce()}, \code{FUN} is
    the callback function to apply to each **viewport** in \code{grid}.
    It must be able to accept as input any of the viewports in \code{grid}.

    For \code{blockReduce()}, \code{init <- FUN(block, init)} will be
    performed on each block so \code{FUN} must take at least two arguments
    (typically \code{block} and \code{init} but the names can differ)
    and should normally return a value of the same type as its 2nd
    argument (\code{init}).

    The same applies for \code{viewportReduce()}, except that
    \code{init <- FUN(viewport, init)} will be performed on each
    **viewport**.
  }
  \item{...}{
    Optional arguments to \code{FUN}.
  }
  \item{grid}{
    An \link{ArrayGrid} object that defines the blocks (or viewports) to
    walk on.

    For \code{blockApply()} and \code{blockReduce()} the supplied grid
    must be compatible with the geometry of \code{x}. If not specified,
    an automatic grid is used. By default \code{\link{defaultAutoGrid}(x)}
    is called to generate an automatic grid. The \emph{automatic grid maker}
    can be changed with \code{\link{setAutoGridMaker}()}.
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
    \pkg{BiocParallel} package for more information about parallel back-ends.
  }
  \item{verbose}{
    Whether block processing progress should be displayed or not.
    If set to \code{NA} (the default), verbosity is controlled
    by \code{DelayedArray:::get_verbose_block_processing()}.
    Setting \code{verbose} to \code{TRUE} or \code{FALSE} overrides this.
  }
  \item{init}{
    The value to pass to the first call to \code{FUN(block, init)}
    (or \code{FUN(viewport, init)}) when \code{blockReduce()}
    (or \code{viewportReduce()}) starts the walk. Note that
    \code{blockReduce()} and \code{viewportReduce()} always operate
    sequentially.
  }
  \item{BREAKIF}{
    An optional callback function that detects a break condition.
    Must return \code{TRUE} or \code{FALSE}.
    At each iteration \code{blockReduce()} (and \code{viewportReduce()})
    will call it on the result of \code{init <- FUN(block, init)}
    (on the result of \code{init <- FUN(viewport, init)} for
    \code{viewportReduce()}) and exit the walk if \code{BREAKIF(init)}
    returned \code{TRUE}.
  }
  \item{envir}{
    Do not use (unless you know what you are doing).
  }
  \item{effective_grid, current_block_id}{
    See Details below.
  }
}

\details{
  \code{effectiveGrid()}, \code{currentBlockId()}, and \code{currentViewport()}
  return the "grid context" for the block/viewport being currently processed.
  By "grid context" we mean:
  \itemize{
    \item The \emph{effective grid}, that is, the user-supplied grid
          or \code{defaultAutoGrid(x)} if the user didn't supply any grid.
    \item The \emph{current block id} (a.k.a. block rank).
    \item The \emph{current viewport}, that is, the \link{ArrayViewport}
          object describing the position of the current block w.r.t. the
          effective grid.
  }
  Note that \code{effectiveGrid()}, \code{currentBlockId()}, and
  \code{currentViewport()} can only be called (with no arguments) from
  **within** the callback functions \code{FUN} and/or \code{BREAKIF}
  passed to \code{blockApply()} and family.

  If you need to be able to test/debug your callback function
  as a standalone function, set an arbitrary \emph{effective grid}
  and \emph{current block id} by calling
  \preformatted{    set_grid_context(effective_grid, current_block_id)}
  **right before** calling the callback function.
}

\value{
  For \code{blockApply()} and \code{viewportApply()}, a list with one
  list element per block/viewport visited.

  For \code{blockReduce()} and \code{viewportReduce()}, the result of
  the last call to \code{FUN}.

  For \code{effectiveGrid()}, the grid (\link{ArrayGrid} object) being
  effectively used.

  For \code{currentBlockId()}, the id (a.k.a. rank) of the current block.

  For \code{currentViewport()}, the viewport (\link{ArrayViewport} object)
  of the current block.
}

\seealso{
  \itemize{
    \item \code{\link{defaultAutoGrid}} and family to generate automatic
          grids to use for block processing of array-like objects.

    \item \link{ArrayGrid} and \link{ArrayViewport} objects.

    \item \code{\link{read_block}} and \code{\link{write_block}}.

    \item \code{\link[BiocParallel]{MulticoreParam}},
          \code{\link[BiocParallel]{SnowParam}}, and
          \code{\link[BiocParallel]{bpparam}}, from the \pkg{BiocParallel}
          package.

    \item \link{DelayedArray} objects.
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
blockApply(m, function(block) sum(log(block + 0.5)), grid=m_grid,
           verbose=TRUE)

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
           BPPARAM=BPPARAM, verbose=TRUE)
## Note that blocks can be visited in any order!

## ---------------------------------------------------------------------
## blockReduce()
## ---------------------------------------------------------------------
FUN <- function(block, init) anyNA(block) || init
blockReduce(FUN, m, init=FALSE, grid=m_grid, verbose=TRUE)

m[10, 1] <- NA
blockReduce(FUN, m, init=FALSE, grid=m_grid, verbose=TRUE)

## With early bailout:
blockReduce(FUN, m, init=FALSE, BREAKIF=identity, grid=m_grid,
            verbose=TRUE)

## Note that this is how the anyNA() method for DelayedArray objects is
## implemented.
}
\keyword{methods}
