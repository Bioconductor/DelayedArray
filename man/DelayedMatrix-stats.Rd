\name{DelayedMatrix-stats}

\alias{DelayedMatrix-stats}

\alias{rowSums}
\alias{rowSums,DelayedMatrix-method}
\alias{colSums}
\alias{colSums,DelayedMatrix-method}

\alias{rowMeans}
\alias{rowMeans,DelayedMatrix-method}
\alias{colMeans}
\alias{colMeans,DelayedMatrix-method}

\alias{rowMins}
\alias{rowMins,DelayedMatrix-method}
\alias{colMins}
\alias{colMins,DelayedMatrix-method}

\alias{rowMaxs}
\alias{rowMaxs,DelayedMatrix-method}
\alias{colMaxs}
\alias{colMaxs,DelayedMatrix-method}

\alias{rowRanges}
\alias{rowRanges,DelayedMatrix-method}
\alias{colRanges}
\alias{colRanges,DelayedMatrix-method}

\title{DelayedMatrix row/col summarization}

\description{
  Only a small number of row/col summarization methods are provided by
  the \pkg{DelayedArray} package.

  See the \pkg{DelayedMatrixStats} package for an extensive set of
  row/col summarization methods.
}

\usage{
\S4method{rowSums}{DelayedMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colSums}{DelayedMatrix}(x, na.rm=FALSE, dims=1)

\S4method{rowMeans}{DelayedMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colMeans}{DelayedMatrix}(x, na.rm=FALSE, dims=1)

\S4method{rowMins}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colMins}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)

\S4method{rowMaxs}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colMaxs}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)

\S4method{rowRanges}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colRanges}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
}

\arguments{
  \item{x}{
    A \link{DelayedMatrix} object.
  }
  \item{na.rm}{
    Should missing values (including \code{NaN}) be omitted from the
    calculations?
  }
  \item{dims, rows, cols}{
    These arguments are not supported. Don't use them.
  }
}

\details{
  All these operations are block-processed.
}

\seealso{
  \itemize{
    \item The \pkg{DelayedMatrixStats} package for more row/col summarization
          methods for \link{DelayedMatrix} objects.

    \item \code{\link[base]{rowSums}} in the \pkg{base} package and
          \code{\link[matrixStats]{rowMaxs}} in the \pkg{matrixStats} package
          for row/col summarization of an ordinary matrix.

    \item \link{DelayedMatrix-utils} for other common operations on
          \link{DelayedMatrix} objects.

    \item \link{DelayedMatrix} objects.

    \item \link[base]{matrix} objects in base R.
  }
}

\examples{
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M12 <- rbind(M1, t(M2))        # delayed

## All these operations are block-processed.

rowSums(M12)
colSums(M12)

rowMeans(M12)
colMeans(M12)

rmins <- rowMins(M12)
cmins <- colMins(M12)

rmaxs <- rowMaxs(M12)
cmaxs <- colMaxs(M12)

rranges <- rowRanges(M12)
cranges <- colRanges(M12)

## Sanity checks:
m12 <- rbind(as.matrix(M1), t(as.matrix(M2)))
stopifnot(identical(rowSums(M12), rowSums(m12)))
stopifnot(identical(colSums(M12), colSums(m12)))
stopifnot(identical(rowMeans(M12), rowMeans(m12)))
stopifnot(identical(colMeans(M12), colMeans(m12)))
stopifnot(identical(rmins, rowMins(m12)))
stopifnot(identical(cmins, colMins(m12)))
stopifnot(identical(rmaxs, rowMaxs(m12)))
stopifnot(identical(cmaxs, colMaxs(m12)))
stopifnot(identical(rranges, cbind(rmins, rmaxs, deparse.level=0)))
stopifnot(identical(cranges, cbind(cmins, cmaxs, deparse.level=0)))
}
\keyword{methods}
