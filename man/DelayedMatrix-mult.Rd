\name{DelayedMatrix-mult}

\alias{DelayedMatrix-mult}

\alias{\%*\%}
\alias{\%*\%,DelayedMatrix,ANY-method}
\alias{\%*\%,ANY,DelayedMatrix-method}
\alias{crossprod}
\alias{crossprod,DelayedMatrix,ANY-method}
\alias{crossprod,ANY,DelayedMatrix-method}
\alias{tcrossprod}
\alias{tcrossprod,DelayedMatrix,ANY-method}
\alias{tcrossprod,ANY,DelayedMatrix-method}

\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}
\alias{crossprod,DelayedMatrix,DelayedMatrix-method}
\alias{tcrossprod,DelayedMatrix,DelayedMatrix-method}
\alias{crossprod,DelayedMatrix,missing-method}
\alias{tcrossprod,DelayedMatrix,missing-method}

\title{DelayedMatrix multiplication and cross-product}

\description{
  Like ordinary matrices in base R, \link{DelayedMatrix} objects and
  derivatives can be multiplied with the \code{\%*\%} operator. They also
  support \code{\link[base]{crossprod}()} and \code{\link[base]{tcrossprod}()}.
}

\details{
  Note that matrix multiplication is not delayed: the output matrix is
  realized block by block. The \emph{automatic realization backend}
  controls where realization happens e.g. in memory as an ordinary matrix
  if not set (i.e. set to \code{NULL}), or in an HDF5 file if set to
  \code{"HDF5Array"}.
  See \code{?\link{setAutoRealizationBackend}} for more information about
  realization backends.
}

\value{
  The object returned by matrix multiplication involving at least one
  \link{DelayedMatrix} object will be either:
  \itemize{
    \item An ordinary matrix if the \emph{automatic realization backend}
          is \code{NULL} (the default).

    \item A \link{DelayedMatrix} object if the \emph{automatic realization
          backend} is not \code{NULL}. In this case, the returned
          \link{DelayedMatrix} object will be either \emph{pristine}
          or made of several \emph{pristine} \link{DelayedMatrix} objects
          bound together (via \code{rbind()} or \code{cbind()}, both
          are delayed operations).

          For example, if the \emph{automatic realization backend} is
          \code{"HDF5Array"}, then the returned \link{DelayedMatrix} object
          will be either an \link[HDF5Array]{HDF5Array} object, or it will
          be a \link{DelayedMatrix} object made of several
          \link[HDF5Array]{HDF5Array} objects bound together.
  }
}

\seealso{
  \itemize{
    \item \code{\link[base]{\%*\%}} and \code{\link[base]{crossprod}} in base R.

    \item \link{DelayedMatrix-stats} for \link{DelayedMatrix} row/col
          summarization.

    \item \link{DelayedMatrix-rowsum} for \code{rowsum()} and
          \code{colsum()} methods for \link{DelayedMatrix} objects.

    \item \code{\link{getAutoRealizationBackend}} and
          \code{\link{setAutoRealizationBackend}} for getting and setting
          the \emph{automatic realization backend}.

    \item \link{DelayedArray} objects.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writing an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.
  }
}

\examples{
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)
M1 <- HDF5Array(toy_h5, "M1")

m <- matrix(runif(50000), ncol=nrow(M1))

## Set backend to NULL for in-memory realization:
setAutoRealizationBackend()
p1 <- m \%*\% M1  # an ordinary matrix

## Set backend to HDF5Array for realization in HDF5 file:
setAutoRealizationBackend("HDF5Array")
P2 <- m \%*\% M1  # an HDF5Array object
P2
path(P2)  # HDF5 file where the result got written

## Sanity checks:
stopifnot(
  is.matrix(p1),
  all.equal(p1, m \%*\% as.matrix(M1)),
  is(P2, "HDF5Array"),
  all.equal(as.matrix(P2), p1)
)
}
\keyword{array}
\keyword{methods}
