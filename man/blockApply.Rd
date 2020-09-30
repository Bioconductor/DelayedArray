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

\title{Block processing of an array-like object}

\description{
  A set of utilities for processing an array-like object block by block.
}

\usage{
blockApply(x, FUN, ..., grid=NULL, as.sparse=FALSE, BPPARAM=getAutoBPPARAM())
blockReduce(FUN, x, init, BREAKIF=NULL, grid=NULL, as.sparse=FALSE)

effectiveGrid(envir=parent.frame(2))
currentBlockId(envir=parent.frame(2))
currentViewport(envir=parent.frame(2))

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
  \item{envir}{
    Do not use (unless you know what you are doing).
  }
  \item{BPPARAM}{
    Coming soon...
  }
  \item{init}{
    Coming soon...
  }
  \item{BREAKIF}{
    Coming soon...
  }
}

\details{
  \code{effectiveGrid()}, \code{currentBlockId()}, and \code{currentViewport()}
  can only be called from within the callback functions passed to
  \code{blockApply()} or \code{blockReduce()}, that is, from within
  functions \code{FUN} and \code{BREAKIF}.
  They can be used to provide the grid context for the block being currently
  processed. They return:
  \itemize{
    \item The grid (\link{ArrayGrid} object) being effectively used for
          \code{effectiveGrid()}.
    \item The id (or rank) of the current block.
    \item The viewport (\link{ArrayViewport} object) of the current block.
  }
}

\seealso{
  \itemize{
    \item \code{\link{defaultAutoGrid}} and family to define grids to use
          in the context of block processing of array-like objects.

    \item \code{\link{read_block}}.

    \item \link{ArrayGrid} and \link{ArrayViewport} objects.

    \item \code{\link[BiocParallel]{MulticoreParam}},
          \code{\link[BiocParallel]{SnowParam}}, and
          \code{\link[BiocParallel]{bpparam}}, from the \pkg{BiocParallel}
          package.

    \item \link{DelayedArray} objects.

    \item Advanced users: \link{RealizationSink} objects for writing an
          array-like object block by block to disk (or to memory).
  }
}

\examples{
## ---------------------------------------------------------------------
## blockApply()
## ---------------------------------------------------------------------

m <- matrix(1:60, ncol=6)
m_grid <- defaultAutoGrid(m, block.length=16, block.shape="hypercube")
blockApply(m, identity, grid=m_grid)
blockApply(m, sum, grid=m_grid)
blockApply(m, function(block) sum(log(block + 0.5)), grid=m_grid)

blockApply(m, function(block) {block + currentBlockId()*1e3}, grid=m_grid)
blockApply(m, function(block) currentViewport(), grid=m_grid)
blockApply(m, dim, grid=m_grid)

## More examples coming soon...

## ---------------------------------------------------------------------
## blockReduce()
## ---------------------------------------------------------------------

## Coming soon...
}
\keyword{methods}
