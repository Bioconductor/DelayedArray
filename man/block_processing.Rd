\name{block_processing}

\alias{block_processing}

\alias{blockApply}
\alias{blockReduce}

\alias{effectiveGrid}
\alias{currentBlockId}
\alias{currentViewport}

\alias{read_block}
\alias{read_block,ANY-method}
\alias{write_block}
\alias{write_block,ANY-method}

\alias{write_array_to_sink}

\alias{getDefaultGridMaker}
\alias{setDefaultGridMaker}

\title{Block processing of an array-like object}

\description{
  A set of utilities for processing an array-like object block by block.
}

\usage{
blockApply(x, FUN, ..., grid=NULL, BPREDO=list(), BPPARAM=bpparam())
blockReduce(FUN, x, init, BREAKIF=NULL, grid=NULL)

effectiveGrid(block)
currentBlockId(block)
currentViewport(block)

read_block(x, viewport)
write_block(x, viewport, block)

getDefaultGridMaker()
setDefaultGridMaker(GRIDMAKER="blockGrid")
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{FUN}{
    Coming soon...
  }
  \item{...}{
    Coming soon...
  }
  \item{grid}{
    Coming soon...
  }
  \item{BPREDO}{
    Coming soon...
  }
  \item{BPPARAM}{
    Coming soon...
  }
  \item{init}{
    Coming soon...
  }
  \item{BREAKIF}{
    Coming soon...
  }
  \item{block}{
    Coming soon...
  }
  \item{viewport}{
    Coming soon...
  }
  \item{GRIDMAKER}{
    The function to use by \code{blockApply()} and \code{blockReduce()}
    as \emph{default grid maker}. If no grid is explicitely supplied to
    \code{blockApply()} or \code{blockReduce()}, a grid will be made by
    calling the \emph{default grid maker} function on array-like object
    \code{x}. This call must return an \link{ArrayGrid} object, say
    \code{grid}, that is compatible with \code{x} i.e. such that
    \code{refdim(grid)} is identical to \code{dim(x)}.

    \code{GRIDMAKER} can be specified as a function or as a single string
    naming a function. It can be a user-defined function or a pre-defined
    grid maker like \code{\link{blockGrid}}, \code{\link{rowGrid}}, or
    \code{\link{colGrid}}.

    The default grid maker is set to \code{\link{blockGrid}} at package
    startup and can be reset anytime to this value by calling
    \code{setDefaultGridMaker()} with no argument.
  }
}

\details{
  Coming soon...
}

\seealso{
  \itemize{
    \item \code{\link{blockGrid}} to define grids to use in the context
          of block processing of array-like objects.

    \item \link{DelayedArray} objects.

    \item Advanced users: \link{RealizationSink} objects for writing an
          array-like object block by block to disk (or to memory).
  }
}

\examples{
## ---------------------------------------------------------------------
## blockApply()
## ---------------------------------------------------------------------

## Coming soon...

## ---------------------------------------------------------------------
## blockReduce()
## ---------------------------------------------------------------------

## Coming soon...

## ---------------------------------------------------------------------
## read_block() / write_block()
## ---------------------------------------------------------------------
## 2 utilities for reading/writing blocks from/to an array-like object.

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
for (b in seq_along(grid0)) {
    viewport <- grid0[[b]]
    block <- read_block(m1, viewport)
    block <- b * 1000L + block
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

## ---------------------------------------------------------------------
## CONTROL THE DEFAULT GRID MAKER
## ---------------------------------------------------------------------
getDefaultGridMaker()
setDefaultGridMaker(function(x) colGrid(x, ncol=5))
getDefaultGridMaker()

m <- matrix(runif(600), ncol=12)
blockApply(m, currentViewport)

## Reset default grid maker to factory settings:
setDefaultGridMaker()
}
\keyword{methods}
