\name{write_block}

\alias{write_block}
\alias{write_block,ANY-method}

\alias{class:RealizationSink}
\alias{RealizationSink-class}
\alias{RealizationSink}

\alias{close,RealizationSink-method}

\alias{class:arrayRealizationSink}
\alias{arrayRealizationSink-class}

\alias{dim,arrayRealizationSink-method}
\alias{write_block,arrayRealizationSink-method}
\alias{coerce,arrayRealizationSink,DelayedArray-method}

\alias{AutoRealizationSink}

\alias{supportedRealizationBackends}
\alias{getAutoRealizationBackend}
\alias{setAutoRealizationBackend}
\alias{getRealizationBackend}
\alias{setRealizationBackend}

\alias{sinkApply}

\title{Write blocks of data to a RealizationSink}

\description{
  Use \code{write_block()} to write a block of array data to a
  RealizationSink derivative.

  \code{sinkApply()} is a convenience function for walking on a
  RealizationSink derivative, typically for filling it with data.

  Note that \code{write_block()} is typically used inside the callback
  function passed to \code{sinkApply}.
}

\usage{
write_block(sink, viewport, block)

## Walk on a RealizationSink derivative:
sinkApply(sink, FUN, ..., grid=NULL, verbose=NA)

## Backend-agnostic RealizationSink constructor:
AutoRealizationSink(dim, dimnames=NULL, type="double", as.sparse=FALSE)

## Get/set the "automatic realization backend":
getAutoRealizationBackend()
setAutoRealizationBackend(BACKEND=NULL)
supportedRealizationBackends()
}

\arguments{
  \item{sink}{
    A **writable** array-like object, typically a RealizationSink derivative.
    Some important notes:
    \itemize{
      \item \link{DelayedArray} objects are NEVER writable, even when they
            don't carry delayed operations (e.g. \link[HDF5Array]{HDF5Array}
            objects from the \pkg{HDF5Array} package), and even when they
            don't carry delayed operations **and** have all their data in
            memory (e.g. \link{RleArray} objects). In other words, there are
            NO exceptions.
      \item RealizationSink is a **virtual** class so \code{sink} will always
            be a RealizationSink derivative, that is, an object that belongs
            to a **concrete** subclass of RealizationSink (e.g. an
            \link[HDF5Array]{HDF5RealizationSink} object from the
            \pkg{HDF5Array} package).
      \item RealizationSink derivatives are considered array-like objects
            i.e. they have dimensions and possibly dimnames.
    }
    Although \code{write_block()} and \code{sinkApply()} will typically be
    used on a RealizationSink derivative, they can also be used on an ordinary
    array or other writable in-memory array-like objects like dgCMatrix objects
    from the \pkg{Matrix} package.
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{sink}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(sink)}.
  }
  \item{block}{
    An ordinary (dense) array or \link{SparseArraySeed} object of the
    same dimensions as \code{viewport}.
  }
  \item{FUN}{
    The callback function to apply to each **viewport** of the grid used
    to walk on \code{sink}. \code{sink <- FUN(sink, viewport, ...)}
    will be performed on each viewport, so \code{FUN} must take at least
    two arguments (typically \code{sink} and \code{viewport} but the names
    can differ) and is expected to return its 1st argument (\code{sink}),
    typically modified (when \code{FUN} contains a call to \code{write_block}).
  }
  \item{...}{
    Additional arguments passed to \code{FUN}.
  }
  \item{grid}{
    The grid used for the walk, that is, an \link{ArrayGrid} object that
    defines the viewports to walk on. It must be compatible with the
    geometry of \code{sink}. If not specified, an automatic grid is
    created by calling \code{\link{defaultSinkAutoGrid}(sink)}, and used.
    See \code{?\link{defaultSinkAutoGrid}} for more information.
  }
  \item{verbose}{
    Whether block processing progress should be displayed or not.
    If set to \code{NA} (the default), verbosity is controlled
    by \code{DelayedArray:::get_verbose_block_processing()}.
    Setting \code{verbose} to \code{TRUE} or \code{FALSE} overrides this.
  }
  \item{dim}{
    The dimensions (specified as an integer vector) of the RealizationSink
    derivative to create.
  }
  \item{dimnames}{
    The dimnames (specified as a list of character vectors or NULLs) of
    the RealizationSink derivative to create.
  }
  \item{type}{
    The type of the data that will be written to the RealizationSink
    derivative to create.
  }
  \item{as.sparse}{
    Whether the data should be written as sparse or not to the
    RealizationSink derivative to create. Not all \emph{realization
    backends} support this.
  }
  \item{BACKEND}{
    \code{NULL} (the default), or a single string specifying the name of
    a realization backend e.g. \code{"HDF5Array"} or \code{"RleArray"}
    etc...
  }
}

\details{
  *** The RealizationSink API ***

  The DelayedArray package provides a simple API for writing blocks
  of array data to disk (or to memory): the "RealizationSink API".
  This API allows the developper to write code that is agnostic about
  the particular on-disk (or in-memory) format being used to store
  the data.

  Here is how to use it:
  \enumerate{
    \item Create a realization sink.
    \item Write blocks of array data to the realization sink with
          one or several calls to \code{write_block()}.
    \item Close the realization sink with \code{close()}.
    \item Coerce the realization sink to \link{DelayedArray}.
  }

  A realization sink is formally represented by a RealizationSink derivative.
  Note that RealizationSink is a virtual class with various concrete
  subclasses like \link[HDF5Array]{HDF5RealizationSink} from the
  \pkg{HDF5Array} package, or \link{RleRealizationSink}.
  Each subclass implements the "RealizationSink API" for a specific
  realization backend.

  To create a realization sink, use the specific constructor function.
  This function should be named as the class itself e.g.
  \code{\link[HDF5Array]{HDF5RealizationSink}()}.

  To create a realization sink in a backend-agnostic way, use
  \code{AutoRealizationSink()}. It will create a RealizationSink derivative
  for the current \emph{automatic realization backend} (see below).

  Once writing to the realization sink is completed, the RealizationSink
  derivative must be closed (with \code{close(sink)}), then coerced to
  \link{DelayedArray} (with \code{as(sink, "DelayedArray")}. What
  specific \link{DelayedArray} derivative this coercion will return
  depends on the specific class of the RealizationSink derivative. For
  example, if \code{sink} is an \link[HDF5Array]{HDF5RealizationSink}
  object from the \pkg{HDF5Array} package, then \code{as(sink, "DelayedArray")}
  will return an \link[HDF5Array]{HDF5Array} instance (the
  \link[HDF5Array]{HDF5Array} class is a \link{DelayedArray} subclass).

  *** The \emph{automatic realization backend} ***

  The \emph{automatic realization backend} is a user-controlled global
  setting that indicates what specific RealizationSink derivative
  \code{AutoRealizationSink()} should return.
  In the context of block processing of a \link{DelayedArray} object,
  this controls where/how realization happens e.g. as an ordinary array
  if not set (i.e. set to \code{NULL}), or as an \link[HDF5Array]{HDF5Array}
  object if set to \code{"HDF5Array"}, or as an \link{RleArray} object
  if set to \code{"RleArray"}, etc...

  Use \code{getAutoRealizationBackend()} or \code{setAutoRealizationBackend()}
  to get or set the \emph{automatic realization backend}.

  Use \code{supportedRealizationBackends()} to get the list of realization
  backends that are currently supported.

  *** Cross realization backend compatibility ***

  Two important things to keep in mind for developers aiming at writing
  code that is compatible across realization backends:
  \itemize{
    \item Realization backends don't necessarily support concurrent
          writing.

          More precisely: Even though it is safe to assume that any
          \link{DelayedArray} object will support concurrent
          \code{read_block()} calls, it is not so safe to assume that
          any RealizationSink derivative will support concurrent calls
          to \code{write_block()}. For example, at the moment,
          \link[HDF5Array]{HDF5RealizationSink} objects do not
          support concurrent writing.

          This means that in order to remain compatible across realization
          backends, code that contains calls to \code{write_block()} should
          NOT be parallelized.

    \item Some realization backends are "linear write only", that is,
          they don't support \emph{random write access}, only
          \emph{linear write access}.

          Such backends will provide a relization sink where the blocks
          of data must be written in linear order (i.e. by ascending rank).
          Furthermore, the geometry of the blocks must also be
          compatible with \emph{linear write access}, that is, they must
          have a "first-dim-grows-first" shape. Concretely this means
          that the grid used to walk on the relization sink must be
          created with something like:
          \preformatted{    colAutoGrid(sink)}
          for a two-dimensional sink, or with something like:
          \preformatted{    defaultSinkAutoGrid(sink)}
          for a sink with an arbitrary number of dimensions.

          See \code{?\link{defaultSinkAutoGrid}} for more information.

          For obvious reasons, "linear write only" realization
          backends do not support concurrent writing.
  }
}

\value{
  For \code{write_block()}, the modified array-like object \code{sink}.

  For \code{sinkApply()}, its 1st argument (\code{sink}), typically modified
  (when callback function \code{FUN} contains a call to \code{write_block}).

  For \code{AutoRealizationSink()}, a RealizationSink derivative with the
  class associated with the current \emph{automatic realization backend}.

  For \code{getAutoRealizationBackend}, \code{NULL} (no backend set yet)
  or a single string specifying the name of the \emph{automatic realization
  backend} currently in use.

  For \code{supportedRealizationBackends}, a data frame with 1 row
  per supported realization backend.
}

\seealso{
  \itemize{
    \item \code{\link{defaultSinkAutoGrid}} to create an automatic grid
          on a RealizationSink derivative.

    \item \link{ArrayGrid} for the formal representation of grids and
          viewports.

    \item \code{\link{read_block}}.

    \item \link{SparseArraySeed} objects.

    \item \code{\link{blockApply}} and family for convenient block
          processing of an array-like object.

    \item \link[HDF5Array]{HDF5RealizationSink} objects in the
          \pkg{HDF5Array} package.

    \item \link[HDF5Array]{HDF5-dump-management} in the \pkg{HDF5Array}
          package to control the location and physical properties of
          automatically created HDF5 datasets.

    \item \link{RleArray} objects.

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## USING THE "RealizationSink API": EXAMPLE 1
## ---------------------------------------------------------------------

## -- STEP 1 --
## Create a realization sink. Note that instead of creating a
## realization sink by calling a backend-specific sink constructor
## (e.g. HDF5Array::HDF5RealizationSink), we set the "automatic
## realization backend" to "HDF5Array" and use backend-agnostic
## constructor AutoRealizationSink():
setAutoRealizationBackend("HDF5Array")
sink <- AutoRealizationSink(c(35L, 50L, 8L))
dim(sink)

## -- STEP 2 --
## Define the grid of viewports to walk on. Here we define a grid made
## of very small viewports on 'sink'. Note that, for real-world use cases,
## block processing will typically use grids made of much bigger
## viewports, usually obtained with defaultSinkAutoGrid().
## Also please note that this grid would not be compatible with "linear
## write only" realization backends. See "Cross realization backend
## compatibility" above in this man page for more information.
sink_grid <- RegularArrayGrid(dim(sink), spacings=c(20, 20, 4))

## -- STEP 3 --
## Walk on the grid, and, for each of viewport, write random data to it.
for (bid in seq_along(sink_grid)) {
    viewport <- sink_grid[[bid]]
    block <- array(runif(length(viewport)), dim=dim(viewport))
    sink <- write_block(sink, viewport, block)
}

## -- An alternative to STEP 3 --
FUN <- function(sink, viewport) {
    block <- array(runif(length(viewport)), dim=dim(viewport))
    write_block(sink, viewport, block)
}
sink <- sinkApply(sink, FUN, grid=sink_grid, verbose=TRUE)

## -- STEP 4 --
## Close the sink and turn it into a DelayedArray object:
close(sink)
A <- as(sink, "DelayedArray")
A

setAutoRealizationBackend()  # unset automatic realization backend

## ---------------------------------------------------------------------
## USING THE "RealizationSink API": EXAMPLE 2
## ---------------------------------------------------------------------

## Say we have a 3D array and want to collapse its 3rd dimension by
## summing the array elements that are stacked vertically, that is, we
## want to compute the matrix M obtained by doing sum(A[i, j, ]) for all
## valid i and j. This is very easy to do with an ordinary array:
collapse_3rd_dim <- function(a) apply(a, MARGIN=1:2, sum)

## or, in a slightly more efficient way:
collapse_3rd_dim <- function(a) {
    m <- matrix(0, nrow=nrow(a), ncol=ncol(a))
    for (z in seq_len(dim(a)[[3]]))
        m <- m + a[ , , z]
    m
}

## With a toy 3D array:
a <- array(runif(8000), dim=c(25, 40, 8))
dim(collapse_3rd_dim(a))
stopifnot(identical(sum(a), sum(collapse_3rd_dim(a))))  # sanity check

## Now say that A is so big that even M wouldn't fit in memory. This is
## a situation where we'd want to compute M block by block:

## -- STEP 1 --
## Create the 2D realization sink:
setAutoRealizationBackend("HDF5Array")
sink <- AutoRealizationSink(dim(a)[1:2])
dim(sink)

## -- STEP 2 --
## Define two grids: one for 'sink' and one for 'a'. Since we're going
## to walk on the two grids simultaneously, read a block from 'a' and
## write it to 'sink', we need to make sure that we define grids that
## are "aligned". More precisely the two grids must have the same number
## of viewports and the viewports in one must correspond to the viewports
## in the other one:
sink_grid <- colAutoGrid(sink, ncol=10)
a_spacings <- c(dim(sink_grid[[1L]]), dim(a)[[3]])
a_grid <- RegularArrayGrid(dim(a), spacings=a_spacings)
dims(sink_grid)  # dimensions of the individual viewports
dims(a_grid)     # dimensions of the individual viewports

## Here is how to check that the two grids are "aligned":
stopifnot(identical(length(sink_grid), length(a_grid)))
stopifnot(identical(dims(sink_grid), dims(a_grid)[ , -3]))

## -- STEP 3 --
## Walk on the two grids simultaneously:
for (bid in seq_along(sink_grid)) {
    ## Read block from 'a'.
    a_viewport <- a_grid[[bid]]
    block <- read_block(a, a_viewport)
    ## Collapse it.
    block <- collapse_3rd_dim(block)
    ## Write the collapsed block to 'sink'.
    sink_viewport <- sink_grid[[bid]]
    sink <- write_block(sink, sink_viewport, block)
}

## -- An alternative to STEP 3 --
FUN <- function(sink, sink_viewport) {
    ## Read block from 'a'.
    bid <- currentBlockId()
    a_viewport <- a_grid[[bid]]
    block <- read_block(a, a_viewport)
    ## Collapse it.
    block <- collapse_3rd_dim(block)
    ## Write the collapsed block to 'sink'.
    write_block(sink, sink_viewport, block)
}
sink <- sinkApply(sink, FUN, grid=sink_grid, verbose=TRUE)

## -- STEP 4 --
## Close the sink and turn it into a DelayedArray object:
close(sink)
M <- as(sink, "DelayedArray")
M

## Sanity check:
stopifnot(identical(collapse_3rd_dim(a), as.array(M)))

setAutoRealizationBackend()  # unset automatic realization backend

## ---------------------------------------------------------------------
## USING THE "RealizationSink API": AN ADVANCED EXAMPLE
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

## 1. Two grids of blocks, one on R and one on G. The blocks in the
##    two grids must contain the same number of columns. We arbitrarily
##    choose to use blocks of 150 columns:
R_grid <- colAutoGrid(R, ncol=150)
G_grid <- colAutoGrid(G, ncol=150)

## 2. The function pFUN(). It will take 2 blocks as input, 1 from R
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

## 3. A RealizationSink derivative where to write the matrices returned
##    by pFUN() as we go:
setAutoRealizationBackend("HDF5Array")
U_sink <- AutoRealizationSink(c(200L, 1000L))

## 4. Finally, we create a grid on U_sink with viewports that contain
##    the same number of columns as the corresponding blocks in R and G:
U_grid <- colAutoGrid(U_sink, ncol=150)

## Note that the three grids should have the same number of viewports:
stopifnot(length(U_grid) == length(R_grid))
stopifnot(length(U_grid) == length(G_grid))

## 5. Now we can proceed. We use a "for" loop to walk on R and G
##    simultaneously, block by block, apply pFUN(), and write the
##    output of pFUN() to U_sink:
for (bid in seq_along(U_grid)) {
    R_block <- read_block(R, R_grid[[bid]])
    G_block <- read_block(G, G_grid[[bid]])
    U_block <- pFUN(R_block, G_block)
    U_sink <- write_block(U_sink, U_grid[[bid]], U_block)
}

## An alternative to the "for" loop is to use sinkApply():
FUN <- function(U_sink, U_viewport) {
    bid <- currentBlockId()
    R_block <- read_block(R, R_grid[[bid]])
    G_block <- read_block(G, G_grid[[bid]])
    U_block <- pFUN(R_block, G_block)
    write_block(U_sink, U_viewport, U_block)
}
U_sink <- sinkApply(U_sink, FUN, grid=U_grid, verbose=TRUE)

close(U_sink)
U <- as(U_sink, "DelayedArray")
U

setAutoRealizationBackend()  # unset automatic realization backend

## ---------------------------------------------------------------------
## VERY BASIC (BUT ALSO VERY ARTIFICIAL) USAGE OF THE
## read_block()/write_block() COMBO
## ---------------------------------------------------------------------

###### On an ordinary matrix ######
m1 <- matrix(1:30, ncol=5)

## Define a viewport on 'm1':
block1_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m1), IRanges(c(3, 2), width=block1_dim))

## Read/tranform/write:
block1 <- read_block(m1, viewport1)
write_block(m1, viewport1, block1 + 1000L)

## Define another viewport on 'm1':
viewport1b <- ArrayViewport(dim(m1), IRanges(c(1, 3), width=block1_dim))

## Read/tranform/write:
write_block(m1, viewport1b, block1 + 1000L)

## No-op:
m <- write_block(m1, viewport1, read_block(m1, viewport1))
stopifnot(identical(m1, m))

########## On a 3D array ##########
a3 <- array(1:60, 5:3)

## Define a viewport on 'a3':
block3_dim <- c(2, 4, 1)
viewport3 <- ArrayViewport(dim(a3), IRanges(c(1, 1, 3), width=block3_dim))

## Read/tranform/write:
block3 <- read_block(a3, viewport3)
write_block(a3, viewport3, block3 + 1000L)

## Define another viewport on 'a3':
viewport3b <- ArrayViewport(dim(a3), IRanges(c(3, 1, 3), width=block3_dim))

## Read/tranform/write:
write_block(a3, viewport3b, block3 + 1000L)

## No-op:
a <- write_block(a3, viewport3, read_block(a3, viewport3))
stopifnot(identical(a3, a))

## ---------------------------------------------------------------------
## LESS BASIC (BUT STILL VERY ARTIFICIAL) USAGE OF THE
## read_block()/write_block() COMBO
## ---------------------------------------------------------------------

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

## ---------------------------------------------------------------------
## supportedRealizationBackends() AND FAMILY
## ---------------------------------------------------------------------

getAutoRealizationBackend()  # no backend set yet

supportedRealizationBackends()
setAutoRealizationBackend("HDF5Array")
getAutoRealizationBackend()  # backend is set to "HDF5Array"
supportedRealizationBackends()

getHDF5DumpChunkLength()
setHDF5DumpChunkLength(500L)
getHDF5DumpChunkShape()

sink <- AutoRealizationSink(c(120L, 50L))
class(sink)  # HDF5-specific realization sink
dim(sink)
chunkdim(sink)

grid <- defaultSinkAutoGrid(sink, block.length=600)
for (bid in seq_along(grid)) {
    viewport <- grid[[bid]]
    block <- 101 * bid + runif(length(viewport))
    dim(block) <- dim(viewport)
    sink <- write_block(sink, viewport, block)
}

close(sink)
A <- as(sink, "DelayedArray")
A

setAutoRealizationBackend()  # unset automatic realization backend
}
\keyword{utilities}
