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

\alias{\%*\%}
\alias{\%*\%,DelayedMatrix,matrix-method}
\alias{\%*\%,matrix,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}

% DelayedMatrix row/col summarization

\alias{rowSums}
\alias{rowSums,DelayedMatrix-method}
\alias{colSums}
\alias{colSums,DelayedMatrix-method}

\alias{rowMeans}
\alias{rowMeans,DelayedMatrix-method}
\alias{colMeans}
\alias{colMeans,DelayedMatrix-method}

\alias{rowMaxs}
\alias{rowMaxs,DelayedMatrix-method}
\alias{colMaxs}
\alias{colMaxs,DelayedMatrix-method}

\alias{rowMins}
\alias{rowMins,DelayedMatrix-method}
\alias{colMins}
\alias{colMins,DelayedMatrix-method}

\alias{rowRanges}
\alias{rowRanges,DelayedMatrix-method}
\alias{colRanges}
\alias{colRanges,DelayedMatrix-method}

% Internal stuff

\alias{subset_seed_as_array,ConformableSeedCombiner-method}

\title{Common operations on DelayedArray objects}

\description{
  Common operations on \link{DelayedArray} objects.
}

\details{
  The operations currently supported on \link{DelayedArray} objects are:

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
    \item matrix multiplication (\%*\%) of an ordinary matrix by a
          \link{DelayedMatrix} object
    \item matrix row/col summarization [\link{DelayedMatrix} objects only]:
          \code{rowSums}, \code{colSums}, \code{rowMeans}, \code{colMeans},
          \code{rowMaxs}, \code{colMaxs}, \code{rowMins}, \code{colMins},
          \code{rowRanges}, and \code{colRanges}
  }
}

\seealso{
  \itemize{
    \item \code{\link[base]{is.na}}, \code{\link[base]{!}},
          \code{\link[base]{mean}}, \code{\link[base]{apply}},
          and \code{\link[base]{\%*\%}} in the \pkg{base} package for the
          corresponding operations on ordinary arrays or matrices.

    \item \code{\link[base]{rowSums}} in the \pkg{base} package and
          \code{\link[matrixStats]{rowMaxs}} in the \pkg{matrixStats} package
          for row/col summarization of an ordinary matrix.

    \item \code{\link{setRealizeBackend}} for how to set a \emph{realize()
          backend}.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writting an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

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

M3 <- rbind(M1, t(M2))
M3

## ---------------------------------------------------------------------
## MATRIX MULTIPLICATION
## ---------------------------------------------------------------------

## Matrix multiplication is not delayed: the output matrix is realized
## block by block. The current "realize() backend" controls where
## realization happens e.g. in memory if set to NULL or in an HDF5 file
## if set to "HDF5Array". See '?realize' for more information about
## "realize() backends".
## The output matrix is returned as a DelayedMatrix object with no delayed
## operations on it. The exact class of the object depends on the backend
## e.g. it will be HDF5Matrix with "HDF5Array" backend.

m <- matrix(runif(50000), ncol=nrow(M1))

## Set backend to NULL for in-memory realization:
setRealizeBackend()
P1 <- m \%*\% M1
P1

## Set backend to HDF5Array for realization in HDF5 file:
setRealizeBackend("HDF5Array")

## With the HDF5Array backend, the output matrix will be written to an
## automatic location on disk:
getHDF5DumpFile()  # HDF5 file where the output matrix will be written
lsHDF5DumpFile()

P2 <- m \%*\% M1
P2

lsHDF5DumpFile()

## Use setHDF5DumpFile() and setHDF5DumpName() from the HDF5Array package
## to control the location of automatically created HDF5 datasets.

stopifnot(identical(as.array(P1), as.array(P2)))

## ---------------------------------------------------------------------
## MATRIX ROW/COL SUMMARIZATION
## ---------------------------------------------------------------------

rowSums(M1)
colSums(M1)

rowMeans(M1)
colMeans(M1)

rmaxs <- rowMaxs(M1)
cmaxs <- colMaxs(M1)

rmins <- rowMins(M1)
cmins <- colMins(M1)

rranges <- rowRanges(M1)
cranges <- colRanges(M1)

stopifnot(identical(cbind(rmins, rmaxs, deparse.level=0), rranges))
stopifnot(identical(cbind(cmins, cmaxs, deparse.level=0), cranges))
}
\keyword{methods}
