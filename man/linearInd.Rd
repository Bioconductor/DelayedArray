\name{linearInd}

\alias{linearInd}

\title{Converting array indices into linear indices}

\description{
  \code{linearInd} performs the reverse conversion of
  \code{base::\link[base]{arrayInd}}, that is, it converts so-called
  \emph{array indices} (i.e. n-uplets) into \emph{linear indices}.
}

\usage{
linearInd(aind, dim)
}

\arguments{
  \item{aind}{
    Typically a numeric matrix like one returned by
    \code{base::\link[base]{arrayInd}}, that is, a matrix where
    each row is an n-uplet representing an array index.
    Each array index must describe a position relative to the
    implicit array i.e. to the array whose dimensions are specified
    via the \code{dim} argument.

    For convenience, \code{aind} can also be specified as a vector with
    one element per dimension in the implicit array, in which case it will
    be treated like a 1-row matrix.

    Note that no bounds checking is performed, that is, values in the j-th
    column of \code{aind} can be < 1 or > \code{dim[j]}.
  }
  \item{dim}{
    An integer vector containing the dimensions of the underlying array.

    Note that \code{dim} can also be an integer matrix, in which case it
    must have the same shape as \code{aind}, that is, 1 row per row in
    \code{aind} and 1 column per dimension.
  }
}

\value{
  An integer vector with one element per row in \code{aind} if
  \code{aind} is a matrix.

  A single integer if \code{aind} is a vector.
}

\seealso{
  \code{\link[base]{arrayInd}} in the \pkg{base} package for the reverse
  conversion.
}

\examples{
dim <- 4:2
linearInd(c(4, 3, 1), dim)
linearInd(c(4, 3, 2), dim)

aind <- rbind(c(1, 1, 1),
              c(2, 1, 1),
              c(3, 1, 1),
              c(4, 1, 1),
              c(1, 2, 1),
              c(1, 1, 2),
              c(4, 3, 2))

linearInd(aind, dim)

## With a matrix of dimensions:

dims <- rbind(c(4L, 3L),
              c(5L, 3L),
              c(6L, 3L))

aind <- rbind(c(1,  2),
              c(1,  2),
              c(1,  2))

linearInd(aind, dims)

## Sanity checks:

dim <- c(33:30, 45L, 30L)
stopifnot(linearInd(rep(1, 6), dim) == 1)
stopifnot(linearInd(dim, dim) == prod(dim))

stopifnot(identical(linearInd(arrayInd(1:120, 6:4), 6:4), 1:120))
stopifnot(identical(linearInd(arrayInd(840:1, 4:7), 4:7), 840:1))
}
\keyword{utilities}
