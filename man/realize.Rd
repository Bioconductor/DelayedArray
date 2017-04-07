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
## Coming soon...
}

\keyword{methods}
