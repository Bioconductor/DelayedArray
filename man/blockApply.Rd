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

effectiveGrid(block)
currentBlockId(block)
currentViewport(block)

getAutoBPPARAM()
setAutoBPPARAM(BPPARAM=NULL)
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{FUN}{
    Coming soon...
  }
  \item{...}{
    Coming soon...
  }
  \item{grid}{
    Coming soon...
  }
  \item{as.sparse}{
    Passed to the internal calls to \code{read_block}.
    See \code{?\link{read_block}} for the details.
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
  \item{block}{
    Coming soon...
  }
}

\details{
  Coming soon...
}

\seealso{
  \itemize{
    \item \code{\link{AutoGrid}} to define grids to use in the context
          of block processing of array-like objects.

    \item \code{\link[BiocParallel]{MulticoreParam}},
          \code{\link[BiocParallel]{SnowParam}}, and
          \code{\link[BiocParallel]{bpparam}}, from the \pkg{BiocParallel}
          package.

    \item \code{\link{read_block}}.

    \item \link{DelayedArray} objects.

    \item Advanced users: \link{RealizationSink} objects for writing an
          array-like object block by block to disk (or to memory).
  }
}

\examples{
## ---------------------------------------------------------------------
## blockApply()
## ---------------------------------------------------------------------

## Coming soon...

## ---------------------------------------------------------------------
## blockReduce()
## ---------------------------------------------------------------------

## Coming soon...
}
\keyword{methods}
