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
\alias{show,ArrayGrid-method}

\alias{isLinear}
\alias{isLinear,ArrayViewport-method}
\alias{isLinear,ArrayGrid-method}

\title{ArrayGrid and ArrayViewport objects}

\description{
  ArrayGrid and ArrayViewport objects are used internally to support
  block processing of array-like objects.
}

\examples{
## ---------------------------------------------------------------------
## ArrayGrid OBJECTS
## ---------------------------------------------------------------------

## Create a regularly-spaced grid on top of a 3700 x 100 x 33 array:
grid <- RegularArrayGrid(c(3700L, 100L, 33L), c(250L, 100L, 10L))

## Dimensions of the reference array:
refdim(grid)

## Number of grid elements along each dimension of the reference array:
dim(grid)

## Total number of grid elements:
length(grid)

## First element in the grid:
grid[[1L]]            # same as grid[[1L, 1L, 1L]]

## Last element in the grid:
grid[[length(grid)]]  # same as grid[[15L, 1L, 4L]]

## Lengths of the grid elements:
lengths(grid)

stopifnot(prod(refdim(grid)) == sum(lengths(grid)))

## ---------------------------------------------------------------------
## ArrayViewport OBJECTS
## ---------------------------------------------------------------------

## Grid elements are ArrayViewport objects:
class(grid[[1L]])

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
grid <- RegularArrayGrid(dim(m0), spacings=c(3L, 2L))
grid
extract_block(m0, grid[[3L]])
extract_block(m0, grid[[1L, 3L]])

## Walk on the grid, colum by column:
m1 <- m0
for (b in seq_along(grid)) {
    viewport <- grid[[b]]
    block <- extract_block(m1, viewport)
    block <- b * 1000L + block
    m1 <- replace_block(m1, viewport, block)
}
m1

## Walk on the grid, row by row:
m2 <- m0
for (i in seq_len(dim(grid)[[1]])) {
  for (j in seq_len(dim(grid)[[2]])) {
    viewport <- grid[[i, j]]
    block <- extract_block(m2, viewport)
    block <- (i * 10L + j) * 1000L + block
    m2 <- replace_block(m2, viewport, block)
  }
}
m2
}
\keyword{internal}
