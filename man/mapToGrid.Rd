\name{mapToGrid}

\alias{mapToGrid}
\alias{mapToGrid,ArbitraryArrayGrid-method}
\alias{mapToGrid,RegularArrayGrid-method}

\alias{mapToRef}
\alias{mapToRef,ArbitraryArrayGrid-method}
\alias{mapToRef,RegularArrayGrid-method}

\title{Map reference array positions to grid positions and vice-versa}

\description{
  Use \code{mapToGrid()} to map a set of reference array positions to
  grid positions.
  Use \code{mapToRef()} for the reverse mapping.
}

\usage{
mapToGrid(aind, grid, linear=FALSE)

mapToRef(major, minor, grid, linear=FALSE)
}

\arguments{
  \item{aind}{
    Typically a numeric matrix like one returned by
    \code{base::\link[base]{arrayInd}}, that is, a matrix where
    each row is an n-uplet representing an \emph{array index}.
    Each array index must describe a position relative to the
    reference array of \code{grid}.

    For convenience, \code{aind} can also be specified as a vector with
    one element per dimension in \code{grid}, in which case it will be
    treated like a 1-row matrix.

    Note that no bounds checking is performed, that is, values in the j-th
    column of \code{aind} can be < 1 or > \code{refdim(grid)[j]}. What
    those values will be mapped to is undefined.
  }
  \item{grid}{
    An ArrayGrid object.
  }
  \item{linear}{
    \code{TRUE} or \code{FALSE}. Controls the format of the output for
    \code{mapToGrid} and the input for \code{mapToRef}.

    By default (i.e. when \code{linear} is \code{FALSE}), the major and minor
    indices returned by \code{mapToGrid} (or taken by \code{mapToRef}) are
    \emph{array indices} (i.e. n-uplets). When \code{linear} is set to
    \code{TRUE}, they are returned (or taken) as \emph{linear indices}.
  }
  \item{major, minor}{
    The \code{major} and \code{minor} components as returned by
    \code{mapToGrid}.
  }
}

\value{
  \itemize{
    \item For \code{mapToGrid()}: A list with 2 components, \code{major}
          and \code{minor}.

          Each row in input matrix \code{aind} is an n-uplet that contains
          the coordinates of a position relative to the reference array of
          \code{grid}.

          By default (i.e. when \code{linear} is \code{FALSE}), the 2
          components of the returned list are integer matrices of the same
          dimensions as the input matrix. A row in the \code{major} (or
          \code{minor}) matrix is called a "major n-uplet" (or "minor n-uplet").
          So for each "input position" (i.e. for each row in the input matrix),
          2 n-uplets are returned: the "major n-uplet" and the "minor n-uplet".
          The "major n-uplet" contains the coordinates of the "input position"
          \emph{in the grid coordinate system}, that is, the coordinates of the
          grid element where the position falls in.
          The "minor n-uplet" represents where exactly the "input position"
          falls \emph{inside} the grid element reported by the "major n-uplet".
          The coordinates in the "minor n-uplet" are \emph{relative} to this
          grid element.

          When \code{linear} is \code{TRUE}, the \code{major} and \code{minor}
          components are returned as linear indices. In this case, both are
          integer vectors containing 1 linear index per "input position".

    \item For \code{mapToRef()}: A numeric matrix like one returned
          by \code{base::\link[base]{arrayInd}} describing positions
          relative to the reference array of \code{grid}.
  }
}

\seealso{
  \itemize{
    \item \link{ArrayGrid} objects.

    \item \code{\link{Lindex2Mindex}} and \code{\link{Mindex2Lindex}} in
          this package (\pkg{DelayedArray}) for converting back and forth
          between \emph{linear indices} and \emph{array indices}.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## Create an arbitrary-spaced grid on top of a 15 x 9 matrix:
grid2 <- ArbitraryArrayGrid(list(c(2L, 7:10, 13L, 15L), c(5:6, 6L, 9L)))

## Create a set of reference array positions:
aind <- rbind(c( 2, 5),  # bottom right corner of 1st grid element
              c( 3, 1),  # top left corner of 2nd grid element
              c(14, 9),  # top right corner of last grid element
              c(15, 7),  # bottom left corner of last grid element
              c(15, 9))  # bottom right corner of last grid element

## Map them to grid positions:
majmin <- mapToGrid(aind, grid2)
majmin

## Reverse mapping:
aind2 <- mapToRef(majmin$major, majmin$minor, grid2)
stopifnot(all.equal(aind2, aind))

majmin <- mapToGrid(aind, grid2, linear=TRUE)
majmin
aind2 <- mapToRef(majmin$major, majmin$minor, grid2, linear=TRUE)
stopifnot(all.equal(aind2, aind))

## Map all the valid positions:
all_positions <- seq_len(prod(refdim(grid2)))
aind <- arrayInd(all_positions, refdim(grid2))
majmin <- data.frame(mapToGrid(aind, grid2, linear=TRUE))
majmin

## Sanity checks:
min_by_maj <- split(majmin$minor,
                    factor(majmin$major, levels=seq_along(grid2)))
stopifnot(identical(lengths(min_by_maj, use.names=FALSE), lengths(grid2)))
stopifnot(all(mapply(isSequence, min_by_maj, lengths(min_by_maj))))
aind2 <- mapToRef(majmin$major, majmin$minor, grid2, linear=TRUE)
stopifnot(identical(aind2, aind))

## More mapping:
grid4 <- RegularArrayGrid(c(50, 20), spacings=c(15L, 9L))
aind <- rbind(c( 1,  1),
              c( 2,  1),
              c( 3,  1),
              c(16,  1),
              c(16,  2),
              c(16, 10),
              c(27, 18))

mapToGrid(aind, grid4)
mapToGrid(aind, grid4, linear=TRUE)
}
\keyword{internal}
