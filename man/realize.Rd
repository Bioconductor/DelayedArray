\name{realize}

\alias{class:RealizationSink}
\alias{RealizationSink-class}

\alias{write_block_to_sink}

\alias{close,RealizationSink-method}

\alias{class:arrayRealizationSink}
\alias{arrayRealizationSink-class}

\alias{dim,arrayRealizationSink-method}
\alias{write_block_to_sink,arrayRealizationSink-method}
\alias{coerce,arrayRealizationSink,DelayedArray-method}

\alias{supportedRealizationBackends}
\alias{getRealizationBackend}
\alias{setRealizationBackend}

\alias{realize}
\alias{realize,ANY-method}

\title{Realize a DelayedArray object}

\description{
  Realize a \link{DelayedArray} object in memory or on disk.
  Get or set the \emph{realization backend} for the current session with
  \code{getRealizationBackend} or \code{setRealizationBackend}.
}

\usage{
supportedRealizationBackends()
getRealizationBackend()
setRealizationBackend(BACKEND=NULL)

realize(x, ...)

\S4method{realize}{ANY}(x, BACKEND=getRealizationBackend())
}

\arguments{
  \item{x}{
    The array-like object to realize.
  }
  \item{...}{
    Additional arguments passed to methods.
  }
  \item{BACKEND}{
    \code{NULL} (the default), or a single string specifying the name of
    the backend. When the backend is set to \code{NULL}, \code{x} is
    realized in memory as an ordinary array by just calling \code{as.array}
    on it.
  }
}

\details{
  The \emph{realization backend} controls where/how realization happens e.g.
  as an ordinary array if set to \code{NULL}, as an \link{RleArray} object
  if set to \code{"RleArray"}, or in an HDF5 file if set to \code{"HDF5Array"}.
}

\value{
  \code{realize(x)} returns a \link{DelayedArray} object. More precisely,
  it returns \code{DelayedArray(as.array(x))} when the backend is set to
  \code{NULL} (the default). Otherwise it returns an instance of the class
  associated with the specified backend (which should extend
  \link{DelayedArray}).
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{RleArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)
M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")
M3 <- rbind(log(M1), t(M2))

supportedRealizationBackends()
getRealizationBackend()  # backend is set to NULL
realize(M3)  # realization as ordinary array

setRealizationBackend("RleArray")
getRealizationBackend()  # backend is set to "RleArray"
realize(M3)  # realization as RleArray object

setRealizationBackend("HDF5Array")
getRealizationBackend()  # backend is set to "HDF5Array"
realize(M3)  # realization in HDF5 file
}

\keyword{methods}
