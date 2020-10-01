\name{read_block}

\alias{read_block}
\alias{read_block,ANY-method}

\alias{is_sparse}
\alias{read_sparse_block}
\alias{read_sparse_block,ANY-method}
\alias{read_sparse_block,SparseArraySeed-method}

\title{Read array blocks}

\description{
  Use \code{read_block} to read a block from an array-like object.
  The function is typically used in the context of block processing
  of array-like objects (typically \link{DelayedArray} objects but
  not necessarily).
}

\usage{
read_block(x, viewport, as.sparse=FALSE)
}

\arguments{
  \item{x}{
    An array-like object, typically a \link{DelayedArray} object
    or derivative.
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{x}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(x)}.
  }
  \item{as.sparse}{
    Can be \code{FALSE}, \code{TRUE}, or \code{NA}.

    If \code{FALSE} (the default), the block is returned as an ordinary
    array (a.k.a. dense array). If \code{TRUE}, it's returned as a
    \link{SparseArraySeed} object. Using \code{as.sparse=NA} is equivalent
    to using \code{as.sparse=is_sparse(x)} and is the most efficient way
    to read a block. (This might become the default in the future.)

    Note that when returned as a 2D \link{SparseArraySeed} object with
    numeric or logical data, a block can easily and efficiently
    be coerced to a \link[Matrix]{sparseMatrix} derivative from the
    \pkg{Matrix} package with \code{as(block, "sparseMatrix")}.
    This will return a dgCMatrix object if \code{type(block)}
    is \code{"double"} or \code{"integer"}, or a lgCMatrix
    object if it's \code{"logical"}.
  }
}

\value{
  The data from \code{x} that belongs to the block delimited by the
  specified viewport. The data is returned as an ordinary (dense) array
  or as a \link{SparseArraySeed} object. In both cases it has the same
  dimensions as the \code{viewport}.
}

\seealso{
  \itemize{
    \item \link{ArrayViewport} objects.

    \item \link{SparseArraySeed} objects.

    \item \code{\link{write_block}}.

    \item \code{\link{blockApply}} for more information about block
          processing of array-like objects.

    \item \code{\link{defaultAutoGrid}} and family to generate automatic
          grids to use for block processing of array-like objects.

    \item \link[Matrix]{dgCMatrix-class} and \link[Matrix]{lgCMatrix-class}
          objects in the \pkg{Matrix} package.

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## TYPICAL USE
## ---------------------------------------------------------------------
## read_block() is typically used in combination with write_block().
## See '?write_block' for typical uses of the read_block/write_block
## combo.

## ---------------------------------------------------------------------
## VERY BASIC (BUT ALSO VERY ARTIFICIAL) EXAMPLE 1:
## Read a block from an ordinary matrix
## ---------------------------------------------------------------------
m1 <- matrix(1:30, ncol=5)
m1

## Define the viewport on 'm1' to read the data from:
block1_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m1), IRanges(c(3, 2), width=block1_dim))
viewport1

## Read the block:
block1 <- read_block(m1, viewport1)  # same as m1[3:6, 2:4, drop=FALSE]
block1

## Sanity checks:
stopifnot(identical(dim(viewport1), dim(block1)))
stopifnot(identical(m1[3:6, 2:4, drop=FALSE], block1))

## ---------------------------------------------------------------------
## VERY BASIC (BUT ALSO VERY ARTIFICIAL) EXAMPLE 2:
## Read a block from a sparse matrix
## ---------------------------------------------------------------------
m2 <- rsparsematrix(12, 20, density=0.2,
                    rand.x=function(n) sample(25, n, replace=TRUE))
m2

## Define the viewport on 'm2' to read the data from:
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

## For comparison, using 'as.sparse=NA' on 'm1' still returns the
## block as an ordinary matrix (a.k.a. dense matrix):
read_block(m1, viewport1, as.sparse=NA)

## Sanity checks:
stopifnot(identical(dim(viewport2), dim(block2)))
stopifnot(identical(dim(viewport2), dim(block2b)))
stopifnot(identical(block2, as.array(block2b)))

## ---------------------------------------------------------------------
## VERY BASIC (BUT ALSO VERY ARTIFICIAL) EXAMPLE 3:
## Read a block from a 3D array
## ---------------------------------------------------------------------
a3 <- array(1:60, 5:3)

## Define the viewport on 'a3' to read the data from:
block3_dim <- c(2, 4, 1)
viewport3 <- ArrayViewport(dim(a3), IRanges(c(1, 1, 3), width=block3_dim))
viewport3

## Read the block:
block3 <- read_block(a3, viewport3)  # same as a3[1:2, 1:4, 3, drop=FALSE]
block3

## Note that unlike [, read_block() never drops dimensions.

## Sanity checks:
stopifnot(identical(dim(viewport3), dim(block3)))
stopifnot(identical(a3[1:2, 1:4, 3, drop=FALSE], block3))
}
\keyword{methods}
