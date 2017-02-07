\name{DelayedArray-utils}

\alias{DelayedArray-utils}

% DelayedArray utils

\alias{dim,ConformableSeedCombiner-method}
\alias{dimnames,ConformableSeedCombiner-method}

\alias{+,DelayedArray,missing-method}
\alias{-,DelayedArray,missing-method}

\alias{pmax2}
\alias{pmin2}
\alias{pmax2,ANY,ANY-method}
\alias{pmin2,ANY,ANY-method}
\alias{pmax2,DelayedArray,vector-method}
\alias{pmin2,DelayedArray,vector-method}
\alias{pmax2,vector,DelayedArray-method}
\alias{pmin2,vector,DelayedArray-method}
\alias{pmax2,DelayedArray,DelayedArray-method}
\alias{pmin2,DelayedArray,DelayedArray-method}

\alias{is.na,DelayedArray-method}
\alias{is.finite,DelayedArray-method}
\alias{is.infinite,DelayedArray-method}
\alias{is.nan,DelayedArray-method}
\alias{!,DelayedArray-method}

\alias{nchar,DelayedArray-method}
\alias{tolower,DelayedArray-method}
\alias{toupper,DelayedArray-method}

\alias{round,DelayedArray-method}
\alias{signif,DelayedArray-method}

\alias{anyNA,DelayedArray-method}
\alias{which,DelayedArray-method}

\alias{mean.DelayedArray}
\alias{mean,DelayedArray-method}

\alias{apply}
\alias{apply,DelayedArray-method}

% DelayedMatrix utils

\alias{rowSums}
\alias{rowSums,DelayedMatrix-method}
\alias{colSums}
\alias{colSums,DelayedMatrix-method}
\alias{rowMeans}
\alias{rowMeans,DelayedMatrix-method}
\alias{colMeans}
\alias{colMeans,DelayedMatrix-method}

\alias{\%*\%}
\alias{\%*\%,DelayedMatrix,matrix-method}
\alias{\%*\%,matrix,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}

% Internal stuff

\alias{subset_seed_as_array,ConformableSeedCombiner-method}

\title{Common operations on DelayedArray objects}

\description{
  Common operations on \link{DelayedArray} objects.
}

\details{
  The operations currently supported by \link{DelayedArray} objects are:

  Delayed operations:
  \itemize{
    \item all the members of the \code{\link[methods]{Ops}},
          \code{\link[methods]{Math}}, and \code{\link[methods]{Math2}} groups
    \item \code{!}
    \item \code{is.na}, \code{is.finite}, \code{is.infinite}, \code{is.nan}
    \item \code{nchar}, \code{tolower}, \code{toupper}
    \item \code{pmax2} and \code{pmin2}
    \item \code{rbind} and \code{cbind} (documented in
          \link[DelayedArray]{cbind})
  }

  Block-processed operations:
  \itemize{
    \item \code{anyNA}, \code{which}
    \item all the members of the \code{\link[methods]{Summary}} group
    \item \code{mean}
    \item \code{apply}
    \item \code{rowSums}, \code{colSums}, \code{rowMeans}, and
          \code{colMeans} [\link{DelayedMatrix} objects only]
    \item matrix multiplication (\%*\%) of an ordinary matrix by a
          \link{DelayedMatrix} object
  }
}

\seealso{
  \itemize{
    \item \link[DelayedArray]{cbind} in this package (\pkg{DelayedArray})
          for binding \link{DelayedArray} objects along their rows or columns.

    \item \code{\link{setRealizeBackend}} for how to set a \emph{realize()
          backend}.

    \item \code{\link[HDF5Array]{writeHDF5Dataset}} in the \pkg{HDF5Array}
          package for writting an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \code{\link[base]{is.na}}, \code{\link[base]{!}},
          \code{\link[base]{mean}}, \code{\link[base]{apply}},
          \code{\link[base]{colSums}}, \code{\link[base]{\%*\%}} in the
          \pkg{base} package for the corresponding operations.

    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package
          for the members of the \code{\link[methods]{Ops}},
          \code{\link[methods]{Math}}, and \code{\link[methods]{Math2}} groups.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
range(M1)
M1 >= 0.5 & M1 < 0.75
log(M1)

M2 <- HDF5Array(toy_h5, "M2")
pmax2(M2, 0)

M <- rbind(M1, t(M2))
M
colMeans(M)

## ---------------------------------------------------------------------
## MATRIX MULTIPLICATION
## ---------------------------------------------------------------------

m <- matrix(runif(60), ncol=12)
M <- DelayedArray(matrix(runif(240), nrow=12))

## Matrix multiplication is not delayed. It will realize the result and
## thus one needs to set a "realize() backend" first if none is set yet.
setRealizeBackend("HDF5Array")

## The HDF5Array backend being set, matrix multiplication will write a
## new HDF5 dataset to disk to an automatic location and return an
## HDF5Matrix object that points to this new dataset.

getHDF5DumpFile()  # HDF5 file where the new dataset will be written
lsHDF5DumpFile()

P <- m \%*\% M
P

getHDF5DumpFile()
lsHDF5DumpFile()

## You can use setHDF5DumpFile() and setHDF5DumpName() defined in the
## HDF5Array package to control the location of automatically created
## HDF5 datasets.
}
\keyword{methods}
