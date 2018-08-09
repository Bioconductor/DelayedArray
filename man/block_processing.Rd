\name{block_processing}

\alias{block_processing}

\alias{blockApply}
\alias{blockReduce}

\alias{effectiveGrid}
\alias{currentBlockId}
\alias{currentViewport}

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

    \item \code{\link{read_block}}.

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
