\name{RealizationSink}

\alias{class:RealizationSink}
\alias{RealizationSink-class}

\alias{close,RealizationSink-method}

\alias{class:arrayRealizationSink}
\alias{arrayRealizationSink-class}

\alias{dim,arrayRealizationSink-method}
\alias{write_block,arrayRealizationSink-method}
\alias{coerce,arrayRealizationSink,DelayedArray-method}

\alias{supportedRealizationBackends}
\alias{getRealizationBackend}
\alias{setRealizationBackend}

\alias{RealizationSink}

\title{RealizationSink objects}

\description{
  Get or set the \emph{realization backend} for the current session with
  \code{getRealizationBackend} or \code{setRealizationBackend}.

  Advanced users: Create a RealizationSink object with the backend-agnostic
  \code{RealizationSink()} constructor. Use this object to write an
  array-like object block by block to disk (or to memory).
}

\usage{
supportedRealizationBackends()
getRealizationBackend()
setRealizationBackend(BACKEND=NULL)

RealizationSink(dim, dimnames=NULL, type="double")
}

\arguments{
   \item{BACKEND}{
     \code{NULL} (the default), or a single string specifying the name of
     a \emph{realization backend}.
   }
   \item{dim}{
     The dimensions (specified as an integer vector) of the RealizationSink
     object to create.
   }
   \item{dimnames}{
     The dimnames (specified as a list of character vectors or NULLs) of
     the RealizationSink object to create.
   }
   \item{type}{
     The type of the data that will be written to the RealizationSink
     object to create.
   }
}

\details{
  The \emph{realization backend} controls where/how realization happens e.g.
  as an ordinary array if set to \code{NULL}, as an \link{RleArray} object
  if set to \code{"RleArray"}, or as an \link[HDF5Array]{HDF5Array} object
  if set to \code{"HDF5Array"}.
}

\value{
  \code{supportedRealizationBackends}: A data frame with 1 row per
  supported \emph{realization backend}.

  \code{getRealizationBackend}: \code{NULL} or a single string specifying
  the name of the \emph{realization backend} currently in use.

  \code{RealizationSink}: A RealizationSink object for the current
  \emph{realization backend}.
}

\seealso{
  \itemize{
    \item \code{\link{write_block}}.

    \item \code{\link{blockGrid}} to define grids to use in the context
          of block processing of array-like objects.

    \item \link{DelayedArray} objects.

    \item \link{RleArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[HDF5Array]{HDF5-dump-management} in the \pkg{HDF5Array}
          package to control the location and physical properties of
          automatically created HDF5 datasets.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. supportedRealizationBackends() AND FAMILY
## ---------------------------------------------------------------------
supportedRealizationBackends()
getRealizationBackend()  # backend is set to NULL

setRealizationBackend("HDF5Array")
supportedRealizationBackends()
getRealizationBackend()  # backend is set to "HDF5Array"

## ---------------------------------------------------------------------
## B. A SIMPLE (AND VERY ARTIFICIAL) RealizationSink() EXAMPLE
## ---------------------------------------------------------------------
getHDF5DumpChunkLength()
setHDF5DumpChunkLength(500L)
getHDF5DumpChunkShape()

sink <- RealizationSink(c(120L, 50L))
dim(sink)
chunkdim(sink)

grid <- blockGrid(sink, block.length=600)
for (bid in seq_along(grid)) {
    viewport <- grid[[bid]]
    block <- 101 * bid + runif(length(viewport))
    dim(block) <- dim(viewport)
    write_block(sink, viewport, block)
}

## Always close the RealizationSink object when you are done writing to
## it and before coercing it to DelayedArray:
close(sink)
A <- as(sink, "DelayedArray")
A

## ---------------------------------------------------------------------
## C. AN ADVANCED EXAMPLE OF USER-IMPLEMENTED BLOCK PROCESSING USING
##    colGrid() AND A REALIZATION SINK
## ---------------------------------------------------------------------
## Say we have 2 matrices with the same number of columns. Each column
## represents a biological sample:
library(HDF5Array)
R <- as(matrix(runif(75000), ncol=1000), "HDF5Array")   # 75 rows
G <- as(matrix(runif(250000), ncol=1000), "HDF5Array")  # 250 rows

## Say we want to compute the matrix U obtained by applying the same
## binary functions FUN() to all samples i.e. U is defined as:
##
##   U[ , j] <- FUN(R[ , j], G[ , j]) for 1 <= j <= 1000
##
## Note that FUN() should return a vector of constant length, say 200,
## so U will be a 200x1000 matrix. A naive implementation would be:
##
##   pFUN <- function(r, g) {
##       stopifnot(ncol(r) == ncol(g))  # sanity check
##       sapply(seq_len(ncol(r)), function(j) FUN(r[ , j], g[ , j]))
##   }
##
## But because U is going to be too big to fit in memory, we can't
## just do pFUN(R, G). So we want to compute U block by block and
## write the blocks to disk as we go. The blocks will be made of full
## columns. Also since we need to walk on 2 matrices at the same time
## (R and G), we can't use blockApply() or blockReduce() so we'll use
## a "for" loop.

## Before we get to the "for" loop, we need 4 things:

## 1) Two grids of blocks, one on R and one on G. The blocks in the
##    2 grids must contain the same number of columns. We arbitrarily
##    choose to use blocks of 150 columns:
R_grid <- colGrid(R, ncol=150)
G_grid <- colGrid(G, ncol=150)

## 2) The function pFUN(). It will take 2 blocks as input, 1 from R
##    and 1 from G, apply FUN() to all the samples in the blocks,
##    and return a matrix with one columns per sample:
pFUN <- function(r, g) {
    stopifnot(ncol(r) == ncol(g))  # sanity check
    ## Return a matrix with 200 rows with random values. Completely
    ## artificial sorry. A realistic example would actually need to
    ## apply the same binary function to r[ ,j] and g[ , j] for
    ## 1 <= j <= ncol(r).
    matrix(runif(200 * ncol(r)), nrow=200)
}

## 3) A RealizationSink object where to write the matrices returned
##    by pFUN() as we go. Note that instead of creating a realization
##    sink by calling a backend-specific sink constructor (e.g.
##    HDF5Array:::HDF5RealizationSink), we use the backend-agnostic
##    constructor RealizationSink() and set the current realization
##    backend to HDF5:

setRealizationBackend("HDF5Array")
U_sink <- RealizationSink(c(200L, 1000L))

## 4) Finally, we create a grid on U_sink with blocks that contain the
##    same number of columns as the corresponding blocks in R and G:

U_grid <- colGrid(U_sink, ncol=150)

## Note that the 3 grids should have the same number of blocks:
stopifnot(length(U_grid) == length(R_grid))
stopifnot(length(U_grid) == length(G_grid))

## Now we can procede. We use a "for" loop where we walk on R and G at
## the same time, block by block, apply pFUN(), and write the output
## of pFUN() to U_sink:
for (bid in seq_along(U_grid)) {
    R_block <- read_block(R, R_grid[[bid]])
    G_block <- read_block(G, G_grid[[bid]])
    U_block <- pFUN(R_block, G_block)
    write_block(U_sink, U_grid[[bid]], U_block)
}

close(U_sink)
U <- as(U_sink, "DelayedArray")

## A note about parallelization: even though concurrent block reading
## from the same object is supported, concurrent writing to a sink is
## not supported yet. So the above code cannot be parallelized at the
## moment.
}
\keyword{utilities}

