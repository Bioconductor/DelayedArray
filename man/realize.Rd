\name{realize}

\alias{supportedRealizeBackends}
\alias{getRealizeBackend}
\alias{setRealizeBackend}

\alias{realize}
\alias{realize,ANY-method}

\title{Realize a DelayedArray object}

\description{
  Realize a \link{DelayedArray} object in memory or on disk.
}

\usage{
supportedRealizeBackends()
getRealizeBackend()
setRealizeBackend(BACKEND=NULL)

realize(x, ...)

\S4method{realize}{ANY}(x, BACKEND=getRealizeBackend())
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
  Coming soon...
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

supportedRealizeBackends()
getRealizeBackend()  # backend is set to NULL
m3 <- realize(M3)    # in-memory realization
class(m3)
dim(m3)

setRealizeBackend("HDF5Array")
getRealizeBackend()  # backend is set to "HDF5Array"
M3 <- realize(M3)    # realization in HDF5 file
M3
}

\keyword{methods}
