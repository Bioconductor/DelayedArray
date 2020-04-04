\name{read_block}

\alias{read_block}
\alias{read_block,ANY-method}
\alias{write_block}
\alias{write_block,ANY-method}

\alias{is_sparse}
\alias{read_sparse_block}
\alias{read_sparse_block,ANY-method}
\alias{read_sparse_block,SparseArraySeed-method}
\alias{write_sparse_block}
\alias{write_sparse_block,ANY-method}

\title{Read/write blocks of array data}

\description{
  2 utilities for reading/writing blocks from/to an array-like object.
}

\usage{
read_block(x, viewport)
write_block(x, viewport, block)
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{viewport}{
    An \link{ArrayViewport} object.
  }
  \item{block}{
    An ordinary array of the same dimensions as \code{viewport}.
  }
}

\seealso{
  \itemize{
    \item \link{ArrayViewport} objects.

    \item \code{\link{blockGrid}} to define grids to use in the context
          of block processing of array-like objects.

    \item \code{\link{blockApply}} for more information about block processing
          of array-like objects.

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
m0 <- matrix(1:30, ncol=5)

block_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m0), IRanges(c(3, 2), width=block_dim))
viewport1

block1 <- read_block(m0, viewport1)
block1

## No-op:
write_block(m0, viewport1, block1)
stopifnot(identical(m0, write_block(m0, viewport1, block1)))

write_block(m0, viewport1, block1 + 100L)

viewport2 <- ArrayViewport(dim(m0), IRanges(c(1, 3), width=block_dim))
write_block(m0, viewport2, block1 + 100L)

## Using a grid:
grid0 <- RegularArrayGrid(dim(m0), spacings=c(3L, 2L))
grid0
length(grid0)  # number of blocks defined by the grid
read_block(m0, grid0[[3L]])  # read 3rd block
read_block(m0, grid0[[1L, 3L]])

## Walk on the grid, colum by column:
m1 <- m0
for (bid in seq_along(grid0)) {
    viewport <- grid0[[bid]]
    block <- read_block(m1, viewport)
    block <- bid * 1000L + block
    m1 <- write_block(m1, viewport, block)
}
m1

## Walk on the grid, row by row:
m2 <- m0
for (i in seq_len(dim(grid0)[[1]])) {
  for (j in seq_len(dim(grid0)[[2]])) {
    viewport <- grid0[[i, j]]
    block <- read_block(m2, viewport)
    block <- (i * 10L + j) * 1000L + block
    m2 <- write_block(m2, viewport, block)
  }
}
m2
}
\keyword{methods}
