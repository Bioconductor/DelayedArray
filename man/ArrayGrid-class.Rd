\name{ArrayGrid-class}
\docType{class}

\alias{class:ArrayViewport}
\alias{ArrayViewport-class}
\alias{ArrayViewport}

\alias{refdim}
\alias{refdim,ArrayViewport-method}
\alias{ranges,ArrayViewport-method}
\alias{start,ArrayViewport-method}
\alias{width,ArrayViewport-method}
\alias{end,ArrayViewport-method}
\alias{dim,ArrayViewport-method}

\alias{show,ArrayViewport-method}

\alias{makeNindexFromArrayViewport}

\alias{class:ArrayGrid}
\alias{ArrayGrid-class}
\alias{ArrayGrid}

\alias{class:ArbitraryArrayGrid}
\alias{ArbitraryArrayGrid-class}
\alias{ArbitraryArrayGrid}

\alias{class:RegularArrayGrid}
\alias{RegularArrayGrid-class}
\alias{RegularArrayGrid}

\alias{refdim,ArbitraryArrayGrid-method}
\alias{refdim,RegularArrayGrid-method}
\alias{dim,ArbitraryArrayGrid-method}
\alias{dim,RegularArrayGrid-method}
\alias{as.character.ArrayGrid}
\alias{as.character,ArrayGrid-method}
\alias{dims}
\alias{dims,ArrayGrid-method}
\alias{lengths,ArrayGrid-method}
\alias{maxlength}
\alias{maxlength,ANY-method}
\alias{maxlength,ArbitraryArrayGrid-method}
\alias{maxlength,RegularArrayGrid-method}
\alias{show,ArrayGrid-method}

\alias{mapToGrid}
\alias{mapToGrid,ArbitraryArrayGrid-method}
\alias{mapToGrid,RegularArrayGrid-method}

\alias{downsample}
\alias{downsample,ArbitraryArrayGrid-method}
\alias{downsample,RegularArrayGrid-method}

\alias{isLinear}
\alias{isLinear,ArrayViewport-method}
\alias{isLinear,ArrayGrid-method}

\title{ArrayGrid and ArrayViewport objects}

\description{
  ArrayGrid and ArrayViewport objects are used internally to support
  block processing of an array-like object.
}

\usage{
ArbitraryArrayGrid(tickmarks)
RegularArrayGrid(refdim, spacings=refdim)

mapToGrid(aind, grid, linear=FALSE)

downsample(x, ratio=1L)
}

\arguments{
  \item{tickmarks}{
    A list of integer vectors, one along each dimension of the reference
    array, representing the tickmarks along that dimension. Each integer
    vector must be sorted in ascending order. NAs or negative values are
    not allowed.
  }
  \item{refdim}{
    An integer vector containing the dimensions of the reference array.
  }
  \item{spacings}{
    An integer vector specifying the grid spacing along each dimension.
  }
  \item{aind}{
    Typically a numeric matrix like one returned by
    \code{base::\link[base]{arrayInd}}, that is, a matrix where
    each row is an n-uplet representing an array index.
    Each array index must describe a position relative to the
    reference array of \code{grid}.

    For convenience, \code{aind} can also be specified as a vector with
    one element per dimension in \code{grid}, in which case it will be
    treated like a 1-row matrix.

    Note that no bounds checking is performed, that is, values in the j-th
    column of \code{aind} can be < 1 or > \code{refdim(grid)[j]}. What
    those values will be mapped to is undefined.
  }
  \item{grid, x}{
    An ArrayGrid object.
  }
  \item{linear}{
    \code{TRUE} or \code{FALSE}. By default (i.e. when \code{linear} is
    \code{FALSE}), the major and minor indices are returned as array
    indices (i.e. as n-uplets). When \code{linear} is set to \code{TRUE},
    they are returned as linear indices.
  }
  \item{ratio}{
    An integer vector specifying the ratio of the downsampling along
    each dimension. Can be of length 1, in which case the same ratio is
    used along all the dimensions.
  }
}

\value{
  \itemize{
    \item For \code{ArbitraryArrayGrid()}: An ArbitraryArrayGrid instance.

    \item For \code{RegularArrayGrid()}: A RegularArrayGrid instance.

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

    \item For \code{downsample()}: An ArrayGrid object on the same reference
          array than \code{x}.
  }
}

\seealso{
  \itemize{
    \item \link{block_processing} for more information about block-processing
          of an array-like object.

    \item \code{\link{chunkGrid}}.

    \item \link{DelayedArray} objects.

    \item \code{\link{linearInd}} in this package (\pkg{DelayedArray})
          and \code{\link[base]{arrayInd}} in the \pkg{base} package
          for converting between array indices and linear indices.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. ArrayGrid OBJECTS
## ---------------------------------------------------------------------

## Create a regularly-spaced grid on top of a 3700 x 100 x 33 array:
grid1 <- RegularArrayGrid(c(3700L, 100L, 33L), c(250L, 100L, 10L))

## Dimensions of the reference array:
refdim(grid1)

## Number of grid elements along each dimension of the reference array:
dim(grid1)

## Total number of grid elements:
length(grid1)

## First element in the grid:
grid1[[1L]]             # same as grid1[[1L, 1L, 1L]]

## Last element in the grid:
grid1[[length(grid1)]]  # same as grid1[[15L, 1L, 4L]]

## Dimensions of the grid elements:
dims(grid1)             # one row per grid element

## Lengths of the grid elements:
lengths(grid1)          # same as rowProds(dims(grid1))
stopifnot(sum(lengths(grid1)) == prod(refdim(grid1)))

maxlength(grid1)        # does not need to compute lengths(grid1)) first
                        # so is more efficient than max(lengths(grid1))
stopifnot(maxlength(grid1) == max(lengths(grid1)))

## Create an arbitrary-spaced grid on top of a 15 x 9 matrix:
grid2 <- ArbitraryArrayGrid(list(c(2L, 7:10, 13L, 15L), c(5:6, 6L, 9L)))
refdim(grid2)
dim(grid2)
length(grid2)
grid2[[1L]]             # same as grid2[[1L, 1L]]
grid2[[length(grid2)]]  # same as grid2[[15L, 9L]]
lengths(grid2)
stopifnot(sum(lengths(grid2)) == prod(refdim(grid2)))

maxlength(grid2)        # does not need to compute lengths(grid2)) first
                        # so is more efficient than max(lengths(grid2))
stopifnot(maxlength(grid2) == max(lengths(grid2)))

## ---------------------------------------------------------------------
## B. ArrayViewport OBJECTS
## ---------------------------------------------------------------------

## Grid elements are ArrayViewport objects:
class(grid1[[1L]])

m0 <- matrix(1:30, ncol=5)

block_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m0), IRanges(c(3, 2), width=block_dim))
viewport1

dim(viewport1)     # 'block_dim'
length(viewport1)
ranges(viewport1)

## 2 utilities (not exported yet) for extracting/replacing blocks from/in
## an array-like object:
extract_block <- DelayedArray:::extract_block
replace_block <- DelayedArray:::replace_block

block1 <- extract_block(m0, viewport1)
block1

## No-op:
replace_block(m0, viewport1, block1)
stopifnot(identical(m0, replace_block(m0, viewport1, block1)))

replace_block(m0, viewport1, block1 + 100L)

viewport2 <- ArrayViewport(dim(m0), IRanges(c(1, 3), width=block_dim))
replace_block(m0, viewport2, block1 + 100L)

## Using a grid:
grid3 <- RegularArrayGrid(dim(m0), spacings=c(3L, 2L))
grid3
extract_block(m0, grid3[[3L]])
extract_block(m0, grid3[[1L, 3L]])

## Walk on the grid, colum by column:
m1 <- m0
for (b in seq_along(grid3)) {
    viewport <- grid3[[b]]
    block <- extract_block(m1, viewport)
    block <- b * 1000L + block
    m1 <- replace_block(m1, viewport, block)
}
m1

## Walk on the grid, row by row:
m2 <- m0
for (i in seq_len(dim(grid3)[[1]])) {
  for (j in seq_len(dim(grid3)[[2]])) {
    viewport <- grid3[[i, j]]
    block <- extract_block(m2, viewport)
    block <- (i * 10L + j) * 1000L + block
    m2 <- replace_block(m2, viewport, block)
  }
}
m2

## ---------------------------------------------------------------------
## C. mapToGrid()
## ---------------------------------------------------------------------

## Map positions relative to the reference array of grid2 to grid2:
aind <- rbind(c( 2, 5),  # bottom right corner of 1st grid element
              c( 3, 1),  # top left corner of 2nd grid element
              c(14, 9),  # top right corner of last grid element
              c(15, 7),  # bottom left corner of last grid element
              c(15, 9))  # bottom right corner of last grid element

mapToGrid(aind, grid2)
mapToGrid(aind, grid2, linear=TRUE)

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

## ---------------------------------------------------------------------
## D. DOWNSAMPLING AN ArrayGrid OBJECT
## ---------------------------------------------------------------------
## The elements (ArrayViewport) of an ArrayGrid object can be replaced
## with bigger elements obtained by merging adjacent elements. How many
## adjacent elements to merge along each dimension is specified via the
## 'ratio' vector (one integer per dimension). We call this operation
## "downsampling. It can be seen as reducing the "resolution" of a grid
## by the specified ratio (if we think of the grid elements as pixels).
downsample(grid2, 2)
downsample(grid2, 3)
downsample(grid2, 4)

## Downsampling preserves the dimensions of the reference array:
stopifnot(identical(refdim(downsample(grid2, 2)), refdim(grid2)))
stopifnot(identical(refdim(downsample(grid2, 3)), refdim(grid2)))
stopifnot(identical(refdim(downsample(grid2, 4)), refdim(grid2)))

## A big enough ratio will eventually produce the coarsest possible grid
## i.e. a grid with a single grid element covering the entire reference
## array:
grid5 <- downsample(grid2, 7)
length(grid5)
grid5[[1L]]
stopifnot(identical(dim(grid5[[1L]]), refdim(grid5)))

## Downsampling by a ratio of 1 is a no-op:
stopifnot(identical(downsample(grid2, 1), grid2))

## Using one ratio per dimension:
downsample(grid2, c(2, 1))
}
\keyword{internal}
