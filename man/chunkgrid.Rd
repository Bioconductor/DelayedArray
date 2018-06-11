\name{chunkgrid}

\alias{class:integer_OR_NULL}
\alias{integer_OR_NULL-class}
\alias{integer_OR_NULL}

\alias{chunkdim}
\alias{chunkdim,ANY-method}
\alias{chunkdim,DelayedUnaryOp-method}
\alias{chunkdim,DelayedSubset-method}
\alias{chunkdim,DelayedAperm-method}

\alias{chunkgrid}
\alias{chunkgrid,ANY-method}

\title{chunkgrid}

\description{
  \code{chunkgrid} and \code{chunkdim} are internal generic functions not
  aimed to be used directly by the user.
}

\usage{
chunkgrid(x)
chunkdim(x)
}

\arguments{
  \item{x}{
    An array-like object.
  }
}

\details{
  Coming soon...
}

\value{
  \code{chunkgrid} returns NULL or an \link{ArrayGrid} object defining a
  grid on reference array \code{x}.

  \code{chunkdim} returns NULL or the chunk dimensions in an integer vector
  parallel to \code{dim(x)}.
}

\seealso{
  \itemize{
    \item \link{block_processing} for more information about block-processing
          of an array-like object.

    \item \link{DelayedArray} objects.

    \item \link{ArrayGrid} objects.
  }
}

\examples{
## Coming soon...
}
\keyword{internal}
