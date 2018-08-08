\name{makeCappedVolumeBox}

\alias{capped_volume_boxes}

\alias{makeCappedVolumeBox}

\alias{makeRegularArrayGridOfCappedLengthViewports}

\alias{isLinear}
\alias{isLinear,ArrayViewport-method}
\alias{isLinear,ArrayGrid-method}

\title{Utilities to make capped volume boxes}

\description{
  \code{makeCappedVolumeBox} returns the dimensions of the biggest
  multidimensional box (a.k.a. hyperrectangle) that satisfies 3 constraints:
  (1) its volume is capped, (2) it fits in the \emph{constraining box},
  (3) it has the specified shape.

  \code{makeRegularArrayGridOfCappedLengthViewports} makes a
  \link{RegularArrayGrid} object with grid elements that are capped volume
  boxes with the specified constraints.

  These are low-level utilities used internally to support
  \code{\link{blockGrid}} and family.
}

\usage{
makeCappedVolumeBox(maxvol, maxdim, shape=c("hypercube",
                                            "scale",
                                            "first-dim-grows-first",
                                            "last-dim-grows-first"))

makeRegularArrayGridOfCappedLengthViewports(refdim,
                           viewport_len,
                           viewport_shape=c("hypercube",
                                            "scale",
                                            "first-dim-grows-first",
                                            "last-dim-grows-first"))
}

\arguments{
  \item{maxvol}{
    The maximum volume of the box to return.
  }
  \item{maxdim}{
    The dimensions of the constraining box.
  }
  \item{shape}{
    The shape of the box to return.
  }
  \item{refdim}{
    The dimensions of the reference array of the grid to return.
  }
  \item{viewport_len}{
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
          dimensions are specified via \code{maxdim}.
    \item The shape of the box is as close as possible to the requested shape.
  }

  The supported shapes are:
  \itemize{
    \item \code{hypercube}: The box should be as close as possible to an
          \emph{hypercube} (a.k.a. \emph{n-cube}), that is, the ratio
          between its biggest and smallest dimensions should be as close
          as possible to 1.

    \item \code{scale}: The box should have the same proportions as the
          \emph{constraining box}.

    \item \code{first-dim-grows-first}: The box will be grown along its
          1st dimension first, then along its 2nd dimension, etc...

    \item \code{last-dim-grows-first}: Like \code{first-dim-grows-first}
          but starting along the last dimension.
  }
}

\seealso{
  \itemize{
    \item \code{\link{blockGrid}} to define a grid suitable for
          block processing of an array-like object.

    \item \link{ArrayGrid} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## makeCappedVolumeBox()
## ---------------------------------------------------------------------

maxdim <- c(50, 12)  # dimensions of the "constraining box"

## "hypercube" shape:
makeCappedVolumeBox(40, maxdim)
makeCappedVolumeBox(120, maxdim)
makeCappedVolumeBox(125, maxdim)
makeCappedVolumeBox(200, maxdim)

## "scale" shape:
makeCappedVolumeBox(40, maxdim, shape="scale")
makeCappedVolumeBox(160, maxdim, shape="scale")

## "first-dim-grows-first" and "last-dim-grows-first" shapes:
makeCappedVolumeBox(120, maxdim, shape="first-dim-grows-first")
makeCappedVolumeBox(149, maxdim, shape="first-dim-grows-first")
makeCappedVolumeBox(150, maxdim, shape="first-dim-grows-first")

makeCappedVolumeBox(40, maxdim, shape="last-dim-grows-first")
makeCappedVolumeBox(59, maxdim, shape="last-dim-grows-first")
makeCappedVolumeBox(60, maxdim, shape="last-dim-grows-first")

## ---------------------------------------------------------------------
## makeRegularArrayGridOfCappedLengthViewports()
## ---------------------------------------------------------------------

grid1a <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 40)
grid1a
as.list(grid1a)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid1a))
stopifnot(maxlength(grid1a) <= 40)  # sanity check

grid1b <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 40,
                                            "first-dim-grows-first")
grid1b
as.list(grid1b)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid1b))
stopifnot(maxlength(grid1b) <= 40)  # sanity check

grid2a <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 120)
grid2a
as.list(grid2a)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid2a))
stopifnot(maxlength(grid2a) <= 120)  # sanity check

grid2b <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 120,
                                            "first-dim-grows-first")
grid2b
as.list(grid2b)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid2b))
stopifnot(maxlength(grid2b) <= 120)  # sanity check

grid3a <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 200)
grid3a
as.list(grid3a)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid3a))
stopifnot(maxlength(grid3a) <= 200)  # sanity check

grid3b <- makeRegularArrayGridOfCappedLengthViewports(maxdim, 200,
                                            "first-dim-grows-first")
grid3b
as.list(grid3b)  # turn the grid into a list of ArrayViewport objects
table(lengths(grid3b))
stopifnot(maxlength(grid3b) <= 200)  # sanity check
}
\keyword{utilities}
