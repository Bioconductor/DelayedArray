\name{block_processing}

\alias{block_processing}

\alias{blockApply}
\alias{blockReduce}

\alias{defaultGrid}
\alias{effectiveGrid}
\alias{currentBlockId}
\alias{currentViewport}

\alias{write_array_to_sink}

\title{Block processing of an array-like object}

\description{
  A set of utilities for processing an array-like object block by block.
}

\usage{
blockApply(x, FUN, ..., grid=NULL, BPREDO=list(), BPPARAM=bpparam())
blockReduce(FUN, x, init, BREAKIF=NULL, grid=NULL)

defaultGrid(x, max.block.length=NULL)
effectiveGrid(block)
currentBlockId(block)
currentViewport(block)
}

\arguments{
  \item{x}{
    Coming soon...
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
  \item{BPREDO}{
    Coming soon...
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
  \item{max.block.length}{
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
    \item \link{DelayedArray} objects.

    \item \code{\link{chunkgrid}}.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.
  }
}

\examples{
## Coming soon...
}
\keyword{methods}
