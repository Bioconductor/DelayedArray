\name{realize}

\alias{BLOCK_write_to_sink}

\alias{realize}
\alias{realize,ANY-method}

\title{Realize an object in memory or on disk}

\description{
  \code{realize()} is an S4 generic function.

  The default \code{realize()} method handles the array case. It will
  realize the array-like object (typically a \link{DelayedArray} object)
  in memory or on disk, depending on the \emph{realization backend} specified
  via its \code{BACKEND} argument,
}

\usage{
realize(x, ...)

\S4method{realize}{ANY}(x, BACKEND=getAutoRealizationBackend())
}

\arguments{
  \item{x}{
    An array-like object (typically a \link{DelayedArray} object) for
    the default method.

    Other types of objects can be supported via additional methods.
    For example, the \pkg{SummarizedExperiment} package defines a
    method for \link[SummarizedExperiment]{SummarizedExperiment} objects (see
    \code{?`\link[SummarizedExperiment]{realize,SummarizedExperiment-method}`}).
  }
  \item{...}{
    Additional arguments passed to methods.
  }
  \item{BACKEND}{
    \code{NULL} or a single string specifying the name of a
    \emph{realization backend}. By default, the \emph{automatic
    realization backend} will be used. This is the backend returned by
    \code{\link{getAutoRealizationBackend}()}.
  }
}

\details{
  The default \code{realize()} method realizes an array-like object \code{x}
  in memory if \code{BACKEND} is \code{NULL}, otherwise on disk.

  Note that, when \code{BACKEND} is not \code{NULL}, \code{x} gets realized as
  a "pristine" \link{DelayedArray} object (e.g. an \link[HDF5Array]{HDF5Array}
  object), that is, as a \link{DelayedArray} object that carries no delayed
  operations. This means that, if \code{x} is itself a \link{DelayedArray}
  object, then the returned object is another \link{DelayedArray} object
  semantically equivalent to \code{x} where the delayed operations carried
  by \code{x} have been realized.
}

\value{
  A "pristine" \link{DelayedArray} object if \code{BACKEND} is
  not \code{NULL}.

  Otherwise, an ordinary matrix or array, or a \link{SparseArraySeed}
  object.
}

\seealso{
  \itemize{
    \item \code{\link{getAutoRealizationBackend}} and
          \code{\link{setAutoRealizationBackend}} for getting and setting
          the \emph{automatic realization backend}.

    \item \link{DelayedArray} objects.

    \item \link{RleArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## In-memory realization
## ---------------------------------------------------------------------

a <- array(1:24, dim=4:2)
realize(a, BACKEND=NULL)  # no-op

A <- DelayedArray(a)
realize(log(A), BACKEND=NULL)  # same as 'as.array(log(A))'

## Sanity checks:
stopifnot(identical(realize(a, BACKEND=NULL), a))
stopifnot(identical(realize(log(A), BACKEND=NULL), log(a)))

## ---------------------------------------------------------------------
## On-disk realization
## ---------------------------------------------------------------------

library(HDF5Array)

realize(log(A), BACKEND="HDF5Array")  # same as 'as(log(A), "HDF5Array")'

## ---------------------------------------------------------------------
## Omitting the 'BACKEND' argument
## ---------------------------------------------------------------------

## When 'BACKEND' is not specified, the "automatic realization backend"
## is used. This backend is controlled via setAutoRealizationBackend().

toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)
M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")
M3 <- rbind(log(M1), t(M2)) + 0.5
M3

## Set the "automatic realization backend" to NULL for in-memory
## realization (as ordinary array or SparseArraySeed object):
setAutoRealizationBackend(NULL)
m3 <- realize(M3)  # in-memory realization

registeredRealizationBackends()

setAutoRealizationBackend("RleArray")
realize(M3)  # realization as RleArray object

setAutoRealizationBackend("HDF5Array")
realize(M3)  # on-disk realization (as HDF5Array object)

setAutoRealizationBackend()  # restore default (NULL)
}
\keyword{methods}
