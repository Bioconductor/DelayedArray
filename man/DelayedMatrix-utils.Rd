\name{DelayedMatrix-utils}

\alias{DelayedMatrix-utils}

\alias{t}

\alias{rowsum}
\alias{rowsum.DelayedMatrix}
\alias{rowsum,DelayedMatrix-method}

\alias{\%*\%}
\alias{\%*\%,ANY,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,ANY-method}
\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}

\title{Common operations on DelayedMatrix objects}

\description{
  Common operations on \link{DelayedMatrix} objects.
}

\details{
  In addition to the operations supported on \link{DelayedArray} objects,
  \link{DelayedMatrix} objects support the following operations:

  Delayed operations:
  \itemize{
    \item \code{t}
  }

  Block-processed operations:
  \itemize{
    \item \code{rowsum}
    \item matrix multiplication (\%*\%) of an ordinary matrix by a
          \link{DelayedMatrix} object
    \item matrix row/col summarization (see
          \code{?`\link{DelayedMatrix-stats}`})
  }
}

\seealso{
  \itemize{
    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \link{DelayedArray-stats} for statistical functions on
          \link{DelayedArray} objects.

    \item \link{DelayedMatrix-stats} for \link{DelayedMatrix} row/col
          summarization.

    \item \code{\link{setRealizationBackend}} for how to set a
          \emph{realization backend}.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writing an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## rowsum()
## ---------------------------------------------------------------------
library(HDF5Array)
set.seed(123)
m0 <- matrix(runif(14400000), ncol=2250,
             dimnames=list(NULL, sprintf("C\%04d", 1:2250)))
M0 <- writeHDF5Array(m0, chunkdim=c(200, 250))
dimnames(M0) <- dimnames(m0)
group <- sample(90, nrow(M0), replace=TRUE)

rs0 <- rowsum(M0, group)

## Let's see block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
setAutoBlockSize(2e6)
rs <- rowsum(M0, group)
setAutoBlockSize()

stopifnot(all.equal(rs0, rs))

## ---------------------------------------------------------------------
## MATRIX MULTIPLICATION
## ---------------------------------------------------------------------
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)
M1 <- HDF5Array(toy_h5, "M1")

## Matrix multiplication is not delayed: the output matrix is realized
## block by block. The current "realization backend" controls where
## realization happens e.g. in memory if set to NULL or in an HDF5 file
## if set to "HDF5Array". See '?realize' for more information about
## "realization backends".
## The output matrix is returned as a DelayedMatrix object with no delayed
## operations on it. The exact class of the object depends on the backend
## e.g. it will be HDF5Matrix with "HDF5Array" backend.

m <- matrix(runif(50000), ncol=nrow(M1))

## Set backend to NULL for in-memory realization:
setRealizationBackend()
P1 <- m \%*\% M1
P1

## Set backend to HDF5Array for realization in HDF5 file:
setRealizationBackend("HDF5Array")

## With the HDF5Array backend, the output matrix will be written to an
## automatic location on disk:
getHDF5DumpFile()  # HDF5 file where the output matrix will be written
lsHDF5DumpFile()

P2 <- m \%*\% M1
P2

lsHDF5DumpFile()

## Use setHDF5DumpFile() and setHDF5DumpName() from the HDF5Array package
## to control the location of automatically created HDF5 datasets.

stopifnot(identical(dim(P1), dim(P2)),
          all.equal(as.array(P1), as.array(P2)))
}
\keyword{methods}
