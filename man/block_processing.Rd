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

defaultGrid(x, block.maxlength=NULL, chunk.grid=NULL,
               block.shape=c("hypercube", "linear"))

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
  \item{block.maxlength}{
    Coming soon...
  }
  \item{chunk.grid}{
    Coming soon...
  }
  \item{block.shape}{
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

    \item \code{\link{chunkGrid}}.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.
  }
}

\examples{
m <- matrix(runif(500), ncol=10)

defaultGrid(m, block.maxlength=25)
defaultGrid(m, block.maxlength=10)
defaultGrid(m, block.maxlength=20)
defaultGrid(m, block.maxlength=75)
defaultGrid(m, block.maxlength=100)

defaultGrid(m, block.maxlength=25, block.shape="linear")
defaultGrid(m, block.maxlength=10, block.shape="linear")
defaultGrid(m, block.maxlength=20, block.shape="linear")
defaultGrid(m, block.maxlength=75, block.shape="linear")
defaultGrid(m, block.maxlength=100, block.shape="linear")

## More to come soon...
}
\keyword{methods}
