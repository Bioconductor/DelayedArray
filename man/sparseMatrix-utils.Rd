\name{sparseMatrix-utils}

\alias{sparseMatrix-utils}

\alias{rowsum,ANY-method}
\alias{colsum,ANY-method}
\alias{rowsum.dgCMatrix}
\alias{rowsum,dgCMatrix-method}

\title{sparseMatrix utilities}

\description{
   Some utilities to operate natively on sparseMatrix objects (e.g. dgCMatrix
   and lgCMatrix objects) from the \pkg{Matrix} package. Mostly for internal
   use by the \pkg{DelayedArray} package.
}

\usage{
## rowsum() and colsum() S4 generics:

#rowsum(x, group, reorder=TRUE, ...)
#colsum(x, group, reorder=TRUE, ...)

## Default methods:

\S4method{rowsum}{ANY}(x, group, reorder=TRUE, ...)
\S4method{colsum}{ANY}(x, group, reorder=TRUE, ...)

## rowsum() method for dgCMatrix objects:

\S4method{rowsum}{dgCMatrix}(x, group, reorder=TRUE, ...)
}

\arguments{
  \item{x}{
    A numeric matrix-like object.
  }
  \item{group, reorder, ...}{
    See \code{?base::\link[base]{rowsum}} for a description of
    these arguments.
  }
}

\value{
  See \code{?base::\link[base]{rowsum}} for the value returned
  by the default \code{rowsum} method.

  The default \code{colsum} method returns
  \code{t(rowsum(t(x), group, reorder=reorder, ...))}.
}

\seealso{
  \itemize{
    \item \link{DelayedMatrix-utils} in this package for the \code{rowsum}
          and \code{colsum} methods defined for \link{DelayedMatrix} objects.

    \item \code{base::\link[base]{rowsum}} in the \pkg{base} package for
          the default \code{rowsum} method.

    \item \linkS4class{dgCMatrix} objects in the \pkg{Matrix} package.
  }
}

\examples{
m0 <- rsparsematrix(1e5, 800, density=0.15)  # sparse representation
m <- as.matrix(m0)                           # dense representation
group <- sample(20, nrow(m), replace=TRUE)

## 'rowsum(m0)' is about 4x faster than 'rowsum(m)':
rs0 <- rowsum(m0, group)
rs <- rowsum(m, group)

stopifnot(identical(rs0, rs))
}
\keyword{methods}
