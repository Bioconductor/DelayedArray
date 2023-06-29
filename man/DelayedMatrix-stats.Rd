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

\alias{rowVars}
\alias{rowVars,DelayedMatrix-method}
\alias{colVars}
\alias{colVars,DelayedMatrix-method}

\title{DelayedMatrix row/col summarization}

\description{
  Only a small number of row/col summarization methods are provided by
  the \pkg{DelayedArray} package.

  See the \pkg{DelayedMatrixStats} package for an extensive set of
  row/col summarization methods.
}

\usage{
## N.B.: Showing ONLY the col*() methods (usage of row*() methods is
## the same):

\S4method{colSums}{DelayedMatrix}(x, na.rm=FALSE, dims=1)

\S4method{colMeans}{DelayedMatrix}(x, na.rm=FALSE, dims=1)

\S4method{colMins}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, useNames=TRUE)

\S4method{colMaxs}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, useNames=TRUE)

\S4method{colRanges}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)

\S4method{colVars}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, useNames=TRUE)
}

\arguments{
  \item{x}{
    A \link{DelayedMatrix} object.
  }
  \item{na.rm, useNames, center}{
    See man pages for the corresponding generics in the \pkg{MatrixGenerics}
    package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{rowVars}})
    for a description of these arguments.
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
    \item \link{DelayedMatrix} objects.

    \item The \pkg{DelayedMatrixStats} package for more row/col summarization
          methods for \link{DelayedMatrix} objects.

    \item The man pages for the various generic functions defined
          in the \pkg{MatrixGenerics} package e.g.
          \code{MatrixGenerics::\link[MatrixGenerics]{colVars}} etc...

    \item \link{DelayedMatrix-utils} for other common operations on
          \link{DelayedMatrix} objects.

    \item \link[base]{matrix} objects in base R.
  }
}

\examples{
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M12 <- rbind(M1, t(M2))  # delayed

## All these operations are block-processed.

rsums <- rowSums(M12)
csums <- colSums(M12)

rmeans <- rowMeans(M12)
cmeans <- colMeans(M12)

rmins <- rowMins(M12)
cmins <- colMins(M12)

rmaxs <- rowMaxs(M12)
cmaxs <- colMaxs(M12)

rranges <- rowRanges(M12)
cranges <- colRanges(M12)

rvars <- rowVars(M12, center=rmeans)
cvars <- colVars(M12, center=cmeans)

## Sanity checks:
m12 <- rbind(as.matrix(M1), t(as.matrix(M2)))
stopifnot(
  identical(rsums, rowSums(m12)),
  identical(csums, colSums(m12)),
  identical(rmeans, rowMeans(m12)),
  identical(cmeans, colMeans(m12)),
  identical(rmins, rowMins(m12)),
  identical(cmins, colMins(m12)),
  identical(rmaxs, rowMaxs(m12)),
  identical(cmaxs, colMaxs(m12)),
  identical(rranges, cbind(rmins, rmaxs, deparse.level=0)),
  identical(cranges, cbind(cmins, cmaxs, deparse.level=0)),
  all.equal(rvars, rowVars(m12)),
  all.equal(cvars, colVars(m12))
)
}
\keyword{methods}
