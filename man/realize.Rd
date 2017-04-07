\name{realize}

\alias{supportedRealizeBackends}
\alias{getRealizeBackend}
\alias{setRealizeBackend}

\alias{realize}
\alias{realize,ANY-method}

\title{Realize a DelayedArray object}

\description{
  Realize a DelayedArray object in memory or on disk.
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
  When the backend is set to \code{NULL} (the default), \code{realize(x)}
  returns \code{as.array(x)}. Otherwise it returns an array-like object of
  the class associated with the specified backend.
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
