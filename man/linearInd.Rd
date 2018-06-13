\name{linearInd}

\alias{linearInd}

\title{Turn array indices into linear indices}

\description{
  \code{linearInd} performs the reverse transformation of
  \code{base::\link[base]{arrayInd}} i.e. it transforms array
  indices into linear indices.
}

\usage{
linearInd(aind, dim)
}

\arguments{
  \item{aind}{
    A numeric vector or matrix.

    If a vector, it must have one element per dimension in the underlying
    array. It will be treated as a 1-row matrix.

    If a matrix, it must be like one returned by
    \code{base::\link[base]{arrayInd}} i.e. it must have one row per
    array index and one column per dimension in the underlying array.

    Note that no bounds checking is performed, that is, values in the j-th
    column of \code{aind} can be < 1 or > \code{dim[j]}.
  }
  \item{dim}{
    An integer vector containing the dimensions of the underlying array.
  }
}

\value{
  A single integer if \code{aind} is a vector.

  An integer vector with one element per row in \code{aind} if
  \code{aind} is a matrix.
}

\seealso{
  \code{\link[base]{arrayInd}} in the \pkg{base} package for the reverse
  transformation.
}

\examples{
aind <- rbind(c(1, 1, 1),
              c(2, 1, 1),
              c(3, 1, 1),
              c(4, 1, 1),
              c(1, 2, 1),
              c(1, 1, 2),
              c(4, 3, 2))

linearInd(aind, 4:2)

dim <- sample(50L, 6L, replace=TRUE)
stopifnot(linearInd(rep(1, 6), dim) == 1)
stopifnot(linearInd(dim, dim) == prod(dim))

stopifnot(identical(linearInd(arrayInd(1:120, 6:4), 6:4), 1:120))
stopifnot(identical(linearInd(arrayInd(840:1, 4:7), 4:7), 840:1))
}
\keyword{utilities}
