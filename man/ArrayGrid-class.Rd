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
\alias{length,ArrayViewport-method}

\alias{show,ArrayViewport-method}

\alias{makeNindexFromArrayViewport}

\alias{class:ArrayGrid}
\alias{ArrayGrid-class}
\alias{ArrayGrid}

\alias{class:ArrayArbitraryGrid}
\alias{ArrayArbitraryGrid-class}
\alias{ArrayArbitraryGrid}

\alias{class:ArrayRegularGrid}
\alias{ArrayRegularGrid-class}
\alias{ArrayRegularGrid}

\alias{refdim,ArrayArbitraryGrid-method}
\alias{refdim,ArrayRegularGrid-method}
\alias{dim,ArrayArbitraryGrid-method}
\alias{dim,ArrayRegularGrid-method}
\alias{length,ArrayGrid-method}
\alias{[[,ArrayGrid-method}
\alias{as.character.ArrayGrid}
\alias{as.character,ArrayGrid-method}
\alias{lengths,ArrayGrid-method}
\alias{show,ArrayGrid-method}

\title{ArrayGrid and ArrayViewport objects}

\description{
  ArrayGrid and ArrayViewport objects are used internally to support block
  processing of array-like objects.
}

\examples{
grid <- ArrayRegularGrid(c(2700L, 100L, 33L), c(250L, 100L, 10L))
## Number of grid elements along each dimension of the reference array:
dim(grid)
## Total number of grid elements:
length(grid)
## First element in the grid:
grid[[1L]]            # same as grid[[1L, 1L, 1L]]
## Last element in the grid:
grid[[length(grid)]]  # same as grid[[11L, 1L, 4L]]
}
\keyword{internal}
