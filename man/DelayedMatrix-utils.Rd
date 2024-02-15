\name{DelayedMatrix-utils}

\alias{DelayedMatrix-utils}

\alias{t}

\alias{rowsum}
\alias{rowsum.DelayedMatrix}
\alias{rowsum,DelayedMatrix-method}
\alias{colsum}
\alias{colsum,DelayedMatrix-method}

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
    \item \code{rowsum} and \code{colsum}
    \item matrix multiplication (\%*\%) where at least one of the two
          operands is a \link{DelayedMatrix} object
    \item matrix row/col summarization (see
          \code{?`\link{DelayedMatrix-stats}`})
  }

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
    \item \code{\link[base]{rowsum}} in the \pkg{base} package for
          computing column sums across rows of an ordinary matrix
          for each level of a grouping variable.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \link{DelayedArray-stats} for statistical functions on
          \link{DelayedArray} objects.

    \item \link{DelayedMatrix-stats} for \link{DelayedMatrix} row/col
          summarization.

    \item \code{\link{getAutoRealizationBackend}} and
          \code{\link{setAutoRealizationBackend}} for getting and setting
          the \emph{automatic realization backend}.

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
## rowsum() / colsum()
## ---------------------------------------------------------------------
library(HDF5Array)
set.seed(123)
m0 <- matrix(runif(14400000), ncol=2250,
             dimnames=list(NULL, sprintf("C\%04d", 1:2250)))
M0 <- writeHDF5Array(m0, chunkdim=c(200, 250))
dimnames(M0) <- dimnames(m0)

## --- rowsum() ---

group <- sample(90, nrow(M0), replace=TRUE)  # define groups of rows
rs <- rowsum(M0, group)
rs[1:5, 1:8]
rs2 <- rowsum(M0, group, reorder=FALSE)
rs2[1:5, 1:8]

## Let's see block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
setAutoBlockSize(2e6)
rs3 <- rowsum(M0, group)
setAutoBlockSize()
DelayedArray:::set_verbose_block_processing(FALSE)

## Sanity checks:
stopifnot(all.equal(rowsum(m0, group), rs))
stopifnot(all.equal(rowsum(m0, group, reorder=FALSE), rs2))
stopifnot(all.equal(rs, rs3))

## --- colsum() ---

group <- sample(30, ncol(M0), replace=TRUE)  # define groups of cols
cs <- colsum(M0, group)
cs[1:5, 1:7]
cs2 <- colsum(M0, group, reorder=FALSE)
cs2[1:5, 1:7]

## Sanity checks:
stopifnot(all.equal(colsum(m0, group), cs))
stopifnot(all.equal(cs, t(rowsum(t(m0), group))))
stopifnot(all.equal(cs, t(rowsum(t(M0), group))))
stopifnot(all.equal(colsum(m0, group, reorder=FALSE), cs2))
stopifnot(all.equal(cs2, t(rowsum(t(m0), group, reorder=FALSE))))
stopifnot(all.equal(cs2, t(rowsum(t(M0), group, reorder=FALSE))))

## ---------------------------------------------------------------------
## MATRIX MULTIPLICATION
## ---------------------------------------------------------------------
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
\keyword{methods}
