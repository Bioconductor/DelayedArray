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
\alias{lengths,ArrayGrid-method}
\alias{maxlength}
\alias{maxlength,ANY-method}
\alias{maxlength,ArbitraryArrayGrid-method}
\alias{maxlength,RegularArrayGrid-method}
\alias{show,ArrayGrid-method}

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
  \item{x}{
    An ArrayGrid object.
  }
  \item{ratio}{
    An integer vector specifying the ratio of the downsampling along
    each dimension. Can be of length 1, in which case the same ratio is
    used along all the dimensions.
  }
}

\value{
  An ArbitraryArrayGrid instance for \code{ArbitraryArrayGrid()}.

  A RegularArrayGrid instance for \code{RegularArrayGrid()}.

  An ArrayGrid object on the same reference array than \code{x}
  for \code{downsample()}.
}

\seealso{
  \itemize{
    \item \link{block_processing} for more information about block-processing
          of an array-like object.

    \item \code{\link{chunkgrid}}.

    \item \link{DelayedArray} objects.

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

## Lengths of the grid elements:
lengths(grid1)
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
## C. DOWNSAMPLING AN ArrayGrid OBJECT
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
grid4 <- downsample(grid2, 7)
length(grid4)
grid4[[1L]]
stopifnot(identical(dim(grid4[[1L]]), refdim(grid4)))

## Downsampling by a ratio of 1 is a no-op:
stopifnot(identical(downsample(grid2, 1), grid2))

## Using one ratio per dimension:
downsample(grid2, c(2, 1))
}
\keyword{internal}
