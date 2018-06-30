\name{makeCappedVolumeBox}

\alias{capped_volume_boxes}

\alias{makeCappedVolumeBox}
\alias{makeRegularArrayGridOfCappedLengthViewports}

\title{Utilities to make capped volume boxes}

\description{
  \code{makeCappedVolumeBox} returns the dimensions of the biggest
  multidimensional box (a.k.a. hyperrectangle) that satisfies 3 constraints:
  (1) its volume is capped, (2) it fits in the \emph{constraining box},
  (3) it has a specified shape.

  \code{makeRegularArrayGridOfCappedLengthViewports} makes a
  \link{RegularArrayGrid} object with grid elements that are capped volume
  boxes with the specified constraints.
}

\usage{
makeCappedVolumeBox(maxvol, dim, shape=c("hypercube",
                                         "scale",
                                         "first-dim-grows-first",
                                         "last-dim-grows-first"))

makeRegularArrayGridOfCappedLengthViewports(refdim,
                      viewport_maxlen,
                      viewport_shape=c("hypercube",
                                       "scale",
                                       "first-dim-grows-first",
                                       "last-dim-grows-first"))
}

\arguments{
  \item{maxvol}{
    The maximum volume of the box to return.
  }
  \item{dim}{
    The dimensions of the constraining box.
  }
  \item{shape}{
    The shape of the box to return.
  }
  \item{refdim}{
    The dimensions of the reference array of the grid to return.
  }
  \item{viewport_maxlen}{
    The maximum length of the elements (a.k.a. viewports) of the
    grid to return.
  }
  \item{viewport_shape}{
    The shape of the elements (a.k.a. viewports) of the grid to return.
  }
}

\details{
  \code{makeCappedVolumeBox} returns the dimensions of a box that satisfies
  the following constraints:
  \enumerate{
    \item The volume of the box is as close as possibe to (but no bigger
          than) \code{maxvol}.
    \item The box fits in the \emph{constraining box} i.e. in the box whose
          dimensions are specified via \code{dim}.
    \item The shape of the box is as close as possible to the requested shape.
  }
}

\seealso{
  \itemize{
    \item \link{block_processing} for more information about block processing
          of an array-like object.

    \item \link{ArrayGrid} objects.
  }
}

\examples{
## "hypercube" shape:
makeCappedVolumeBox(40, c(50, 12))
makeCappedVolumeBox(120, c(50, 12))
makeCappedVolumeBox(125, c(50, 12))
makeCappedVolumeBox(200, c(50, 12))

## "scale" shape:
makeCappedVolumeBox(40, c(50, 12), shape="scale")
makeCappedVolumeBox(160, c(50, 12), shape="scale")

## "first-dim-grows-first" and "last-dim-grows-first" shapes:
makeCappedVolumeBox(120, c(50, 12), shape="first-dim-grows-first")
makeCappedVolumeBox(149, c(50, 12), shape="first-dim-grows-first")
makeCappedVolumeBox(150, c(50, 12), shape="first-dim-grows-first")

makeCappedVolumeBox(40, c(50, 12), shape="last-dim-grows-first")
makeCappedVolumeBox(59, c(50, 12), shape="last-dim-grows-first")
makeCappedVolumeBox(60, c(50, 12), shape="last-dim-grows-first")
}
\keyword{methods}
