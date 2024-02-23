\name{DelayedMatrix-rowsum}

\alias{DelayedMatrix-rowsum}

\alias{rowsum}
\alias{rowsum.DelayedMatrix}
\alias{rowsum,DelayedMatrix-method}
\alias{colsum}
\alias{colsum,DelayedMatrix-method}

\title{rowsum() and colsum() on a DelayedMatrix object}

\description{
  Like ordinary matrices in base R, \link{DelayedMatrix} objects and
  derivatives support \code{rowsum()} and \code{colsum()}.
}

\details{
  Note that the \code{rowsum()} and \code{colsum()} operations are not
  delayed: the output matrix is realized block by block.
  The \emph{automatic realization backend} controls where realization
  happens e.g. in memory as an ordinary matrix if not set (i.e. set to
  \code{NULL}), or in an HDF5 file if set to \code{"HDF5Array"}.
  See \code{?\link{setAutoRealizationBackend}} for more information
  about realization backends.
}

\value{
  The object returned by the \code{rowsum()} or \code{colsum()} method
  for \link{DelayedMatrix} objects will be either:
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
    \item \code{\link[base]{rowsum}} in base R.

    \item \code{S4Arrays::\link[S4Arrays]{rowsum}} in the \pkg{S4Arrays}
          package for the \code{rowsum()} and \code{colsum()} S4 generic
          functions.

    \item \code{\link{getAutoRealizationBackend}} and
          \code{\link{setAutoRealizationBackend}} for getting and setting
          the \emph{automatic realization backend}.

    \item \link{DelayedMatrix-stats} for \link{DelayedMatrix} row/col
          summarization.

    \item \link{DelayedMatrix-mult} for \link{DelayedMatrix} multiplication
          and cross-product.

    \item \link{DelayedArray} objects.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writing an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.
  }
}

\examples{
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
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
