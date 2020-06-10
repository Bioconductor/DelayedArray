\name{read_block}

\alias{read_block}
\alias{read_block,ANY-method}

\alias{is_sparse}
\alias{read_sparse_block}
\alias{read_sparse_block,ANY-method}
\alias{read_sparse_block,SparseArraySeed-method}

\alias{write_block}
\alias{write_block,ANY-method}

\title{Read/write blocks of array data}

\description{
  2 utilities for reading/writing blocks from/to an array-like object.
  Their primary purpose is to be used in the context of block processing
  of array-like objects.
}

\usage{
read_block(x, viewport, as.sparse=FALSE)
write_block(x, viewport, block)
}

\arguments{
  \item{x}{
    An array-like object.

    Note that even though \code{write_block()} works when \code{x} is an
    ordinary array or matrix, \code{x} will always be a
    \link{RealizationSink} derivative in the context of block
    processing of \link{DelayedArray} objects.
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{x}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(x)}.
  }
  \item{as.sparse}{
    Can be \code{FALSE}, \code{TRUE}, or \code{NA}.

    If \code{FALSE} (the default), the block is returned as an ordinary
    (dense) array. If \code{TRUE}, it's returned as a \link{SparseArraySeed}
    object. Using \code{as.sparse=NA} is equivalent to using
    \code{as.sparse=is_sparse(x)} and is the most efficient way
    to read a block. (This might become the default in the future.)

    Note that when returned as a 2D \link{SparseArraySeed} object with
    numeric or logical data, a block can easily and efficiently
    be coerced to a \link[Matrix]{sparseMatrix} derivative from the
    \pkg{Matrix} package with \code{as(block, "sparseMatrix")}.
    This will return a dgCMatrix object if \code{type(block)}
    is \code{"double"} or \code{"integer"}, or a lgCMatrix
    object if it's \code{"logical"}.
  }
  \item{block}{
    An ordinary (dense) array or \link{SparseArraySeed} object of the
    same dimensions as \code{viewport}.
  }
}

\value{
  An ordinary (dense) array or \link{SparseArraySeed} object of the
  same dimensions as \code{viewport} for \code{read_block()}.

  The modified array-like object \code{x} for \code{write_block()}.
}

\seealso{
  \itemize{
    \item \link{ArrayViewport} objects.

    \item \link{SparseArraySeed} objects.

    \item \link[Matrix]{dgCMatrix-class} and \link[Matrix]{lgCMatrix-class}
          objects in the \pkg{Matrix} package.

    \item \code{\link{blockGrid}} to define grids to use in the context
          of block processing of array-like objects.

    \item \code{\link{blockApply}} for more information about block
          processing of array-like objects.

    \item Advanced users: \link{RealizationSink} objects for writing an
          array-like object block by block to disk (or to memory).

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## read_block() BASIC USAGE
## ---------------------------------------------------------------------

## == Read a block from an ordinary matrix ==

m1 <- matrix(1:30, ncol=5)
m1
block1_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m1), IRanges(c(3, 2), width=block1_dim))
viewport1

block1 <- read_block(m1, viewport1)
block1

## == Read a block from a sparse matrix ==

m2 <- rsparsematrix(12, 20, density=0.2,
                    rand.x=function(n) sample(25, n, replace=TRUE))
m2
block2_dim <- c(2, 20)
viewport2 <- ArrayViewport(dim(m2), IRanges(c(1, 1), width=block2_dim))
viewport2

## By default, read_block() always returns an ordinary matrix or array:
block2 <- read_block(m2, viewport2)
block2

## It is recommended to use 'as.sparse=NA' rather than 'as.sparse=TRUE'
## or 'as.sparse=FALSE' to let read_block() pick up the optimal
## representation:
block2b <- read_block(m2, viewport2, as.sparse=NA)
class(block2b)  # a SparseArraySeed object
as(block2b, "sparseMatrix")

read_block(m1, viewport1, as.sparse=NA)  # still a dense matrix

## == Read a block from a 3D array ==

a3 <- array(1:60, 5:3)
a3
block3_dim <- c(2, 4, 1)
viewport3 <- ArrayViewport(dim(a3), IRanges(c(1, 1, 3), width=block3_dim))
viewport3
block3 <- read_block(a3, viewport3)
block3

## ---------------------------------------------------------------------
## write_block() BASIC USAGE
## ---------------------------------------------------------------------

write_block(m1, viewport1, block1 + 1000L)

viewport1b <- ArrayViewport(dim(m1), IRanges(c(1, 3), width=block1_dim))
write_block(m1, viewport1b, block1 + 1000L)

write_block(a3, viewport3, block3 + 1000L)

viewport3b <- ArrayViewport(dim(a3), IRanges(c(3, 1, 3), width=block3_dim))
write_block(a3, viewport3b, block3 + 1000L)

## No-ops:
m <- write_block(m1, viewport1, block1)
stopifnot(identical(m1, m))
a <- write_block(a3, viewport3, block3)
stopifnot(identical(a3, a))

## ---------------------------------------------------------------------
## WALKING ON A GRID OF BLOCKS
## ---------------------------------------------------------------------
## In the context of block processing, a grid of blocks is defined on
## the array-like object to process. Then the blocks are processed
## separately in a loop. Here we define a grid made of very small blocks
## on 'm1'. Please note that real-world block processing use grids made
## of much bigger blocks that are obtained with blockGrid().
## See ?blockGrid for more information.

grid1 <- RegularArrayGrid(dim(m1), spacings=c(3L, 2L))
grid1
length(grid1)  # number of blocks defined by the grid
read_block(m1, grid1[[3L]])  # read 3rd block
read_block(m1, grid1[[1L, 3L]])

## Walk on the grid, colum by column:
m1a <- m1
for (bid in seq_along(grid1)) {
    viewport <- grid1[[bid]]
    block <- read_block(m1a, viewport)
    block <- bid * 1000L + block
    m1a <- write_block(m1a, viewport, block)
}
m1a

## Walk on the grid, row by row:
m1b <- m1
for (i in seq_len(dim(grid1)[[1]])) {
  for (j in seq_len(dim(grid1)[[2]])) {
    viewport <- grid1[[i, j]]
    block <- read_block(m1b, viewport)
    block <- (i * 10L + j) * 1000L + block
    m1b <- write_block(m1b, viewport, block)
  }
}
m1b
}
\keyword{methods}
