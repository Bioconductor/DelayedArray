\name{cbind-methods}

\alias{cbind-methods}

\alias{dim,SeedBinder-method}
\alias{dimnames,SeedBinder-method}

\alias{rbind}
\alias{rbind,DelayedMatrix-method}
\alias{rbind,DelayedArray-method}

\alias{cbind}
\alias{cbind,DelayedMatrix-method}
\alias{cbind,DelayedArray-method}

\alias{arbind}
\alias{arbind,DelayedArray-method}

\alias{acbind}
\alias{acbind,DelayedArray-method}

% Internal stuff
\alias{subset_seed_as_array,SeedBinder-method}

\title{Bind DelayedArray objects along their rows or columns}

\description{
  Methods for binding DelayedArray objects along their rows or columns.
}

\details{
  \code{rbind}, \code{cbind}, \code{arbind}, \code{acbind} methods are defined
  for \link{DelayedArray} objects. They perform delayed binding along the rows
  (\code{rbind} and \code{arbind}) or columns (\code{cbind} and \code{acbind})
  of the objects passed to them.
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in the \pkg{base} package for
          rbind/cbind'ing ordinary arrays.

    \item \code{\link[IRanges]{acbind}} in the \pkg{IRanges} package for
          arbind/acbind'ing ordinary arrays.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## rbind/cbind
## ---------------------------------------------------------------------
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M <- rbind(M1, t(M2))
M
colMeans(M)

## ---------------------------------------------------------------------
## arbind/acbind
## ---------------------------------------------------------------------
a1 <- array(1:60, c(3, 5, 4),
            dimnames=list(NULL, paste0("M1y", 1:5), NULL))
a2 <- array(101:240, c(7, 5, 4),
            dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL))
a3 <- array(10001:10100, c(5, 5, 4),
            dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))

A1 <- DelayedArray(a1)
A2 <- DelayedArray(a2)
A3 <- DelayedArray(a3)
A <- arbind(A1, A2, A3)
A

## Sanity check:
stopifnot(identical(arbind(a1, a2, a3), as.array(A)))
}
\keyword{methods}
