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

\title{Write array blocks}

\description{
  Use \code{write_block} to write a block of array data to a RealizationSink
  object. The function is typically used in the context of block processing
  of array-like objects (typically \link{DelayedArray} objects but
  not necessarily).
}

\usage{
write_block(x, viewport, block)

## Backend-agnostic RealizationSink constructor:
AutoRealizationSink(dim, dimnames=NULL, type="double", as.sparse=FALSE)

## Get/set the "automatic realization backend":
getAutoRealizationBackend()
setAutoRealizationBackend(BACKEND=NULL)
supportedRealizationBackends()
}

\arguments{
  \item{x}{
    A **writable** array-like object, typically a RealizationSink derivative.
    Some important notes:
    \itemize{
      \item \link{DelayedArray} objects are NEVER writable, even when they
            don't carry delayed operations (e.g. \link[HDF5Array]{HDF5Array}
            objects from the \pkg{HDF5Array} package), even when they don't
            carry delayed operations AND have all their data in memory (e.g.
            \link{RleArray} objects).
      \item RealizationSink is a virtual class so \code{x} must be an
            instance of a **concrete** RealizationSink subclass (e.g.
            an \link[HDF5Array]{HDF5RealizationSink} object from the
            \pkg{HDF5Array} package).
      \item RealizationSink derivatives are considered array-like objects
            i.e. they have dimensions and possibly dimnames.
    }
    Although \code{write_block()} will typically be used on a RealizationSink
    derivative, it can also be used on an ordinary array or other writable
    in-memory array-like arrays like dgCMatrix objects from the \pkg{Matrix}
    package.
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{x}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(x)}.
  }
  \item{block}{
    An ordinary (dense) array or \link{SparseArraySeed} object of the
    same dimensions as \code{viewport}.
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
  \item{as.sparse}{
    Whether the data should be written as sparse or not to the
    RealizationSink object to create. Not all \emph{realization
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

  A realization sink is formally represented by a RealizationSink object.
  Note that RealizationSink is a virtual class with various concrete
  subclasses like \link[HDF5Array]{HDF5RealizationSink} from the
  \pkg{HDF5Array} package, or \link{RleRealizationSink}.
  Each subclass implements the "RealizationSink API" for a specific
  realization backend.

  To create a realization sink, use the specific constructor function.
  This function should be named as the class itself e.g.
  \code{\link[HDF5Array]{HDF5RealizationSink}()}.

  To create a realization sink in a backend-agnostic way, use
  \code{AutoRealizationSink()}. It will create a RealizationSink object
  for the current \emph{automatic realization backend} (see below).

  Once writing to the realization sink is completed, the RealizationSink
  object must be closed (with \code{close(sink)}), then coerced to
  \link{DelayedArray} (with \code{as(sink, "DelayedArray")}. What
  specific \link{DelayedArray} derivative this coercion will return
  depends on the specific class of the RealizationSink object. For
  example, if \code{sink} is an \link[HDF5Array]{HDF5RealizationSink}
  object from the \pkg{HDF5Array} package, then \code{as(sink, "DelayedArray")}
  will return an \link{HDF5Array} instance (the \link{HDF5Array} class
  is a \link{DelayedArray} subclass).

  *** The \emph{automatic realization backend} ***

  The \emph{automatic realization backend} is a user-controlled global
  setting that indicates what specific RealizationSink object
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
}

\value{
  For \code{write_block()}, the modified array-like object \code{x}.

  For \code{AutoRealizationSink()}, a RealizationSink object for the
  current \emph{automatic realization backend}.

  For \code{getAutoRealizationBackend}, \code{NULL} (no backend set yet)
  or a single string specifying the name of the \emph{automatic realization
  backend} currently in use.

  For \code{supportedRealizationBackends}, a data frame with 1 row
  per supported realization backend.
}

\note{
  IMPORTANT NOTE ABOUT PARALLELIZATION: Even though it is safe to
  assume that any \link{DelayedArray} object will support concurrent
  \code{read_block()} calls, it is not so safe to assume that any
  RealizationSink derivative will support concurrent calls to
  \code{write_block()}. For example, at the moment
  \link[HDF5Array]{HDF5RealizationSink} objects do NOT support
  concurrent writing.

  This means that in order to remain backend-agnostic, code that contains
  calls to \code{write_block()} shoult NOT be parallelized.
}

\seealso{
  \itemize{
    \item \link{ArrayViewport} objects.

    \item \link{SparseArraySeed} objects.

    \item \code{\link{read_block}}.

    \item \code{\link{blockApply}} for more information about block
          processing of array-like objects.

    \item \code{\link{defaultAutoGrid}} and family to generate automatic
          grids to use for block processing of array-like objects.

    \item \link{RleArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[HDF5Array]{HDF5-dump-management} in the \pkg{HDF5Array}
          package to control the location and physical properties of
          automatically created HDF5 datasets.

    \item \link{DelayedArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## USING THE "RealizationSink API": EXAMPLE 1
## ---------------------------------------------------------------------

## 1. Create a realization sink (note that we choose to use an HDF5
##    realization sink but other realization backends would work the
##    same way):
library(HDF5Array)
sink <- HDF5RealizationSink(c(25L, 40L, 8L))
dim(sink)

## 2. Define the grid of blocks to walk on. Here we define a grid made
##    of very small blocks on 'sink'. Please note that real-world block
##    processing typically uses grids made of much bigger blocks obtained
##    with defaultAutoGrid(). See ?defaultAutoGrid for more information.
sink_grid <- RegularArrayGrid(dim(sink), spacings=c(12, 10, 2))

## 3. Walk on the grid, and, for each of block, write random data to it:
for (bid in seq_along(sink_grid)) {
    viewport <- sink_grid[[bid]]
    block <- array(runif(length(viewport)), dim=dim(viewport))
    sink <- write_block(sink, viewport, block)
}

## 4. Close the sink and turn it into a DelayedArray object:
close(sink)
A <- as(sink, "DelayedArray")
A

## Note that, except for step 1, the above code is backend-agnostic.
## To also make step 1 backend-agnostic, replace the call to
## HDF5RealizationSink() with a call to AutoRealizationSink().
## Also note that in order to remain backend-agnostic the above code
## should not be parallelized. See "IMPORTANT NOTE ABOUT PARALLELIZATION"
## above in this man page.

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

## 1. Create the 2D realization sink:
sink <- HDF5RealizationSink(dim(a)[1:2])
dim(sink)

## 2. Define two grids: one for 'sink' and one for 'a'. Since we're going
##    to walk on the two grids simultaneously, read a block from 'a' and
##    write it to 'sink', we need to make sure that we define grids that
##    are "aligned". More precisely they must have the same number of
##    blocks and the blocks in one must correspond to the blocks in the
##    other:
sink_grid <- RegularArrayGrid(dim(sink), spacings=c(12, 10))
a_grid <- RegularArrayGrid(dim(a), spacings=c(12, 10, dim(a)[[3]]))
stopifnot(identical(length(sink_grid), length(a_grid)))
dims(sink_grid)  # dimensions of the individual blocks
dims(a_grid)     # dimensions of the individual blocks
stopifnot(identical(dims(sink_grid), dims(a_grid)[ , -3]))  # sanity check

## 3. Walk on the two grids simultaneously:
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

## 4. Close the sink and turn it into a DelayedArray object:
close(sink)
M <- as(sink, "DelayedArray")
M

## Note that, except for step 1, the above code is backend-agnostic.
## To also make step 1 backend-agnostic, replace the call to
## HDF5RealizationSink() with a call to AutoRealizationSink().
## Also note that in order to remain backend-agnostic the above code
## should not be parallelized. See "IMPORTANT NOTE ABOUT PARALLELIZATION"
## above in this man page.

stopifnot(identical(collapse_3rd_dim(a), as.array(M)))  # sanity check

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
##    2 grids must contain the same number of columns. We arbitrarily
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

## 3. A RealizationSink object where to write the matrices returned
##    by pFUN() as we go. Note that instead of creating a realization
##    sink by calling a backend-specific sink constructor (e.g.
##    HDF5Array::HDF5RealizationSink), we set the "automatic realization
##    backend" to HDF5 and use the backend-agnostic constructor
##    AutoRealizationSink():

setAutoRealizationBackend("HDF5Array")
U_sink <- AutoRealizationSink(c(200L, 1000L))

## 4. Finally, we create a grid on U_sink with blocks that contain the
##    same number of columns as the corresponding blocks in R and G:

U_grid <- colAutoGrid(U_sink, ncol=150)

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
    sink <- write_block(U_sink, U_grid[[bid]], U_block)
}

close(U_sink)
U <- as(U_sink, "DelayedArray")
U

## Note that in order to remain backend-agnostic the above code should
## not be parallelized. See "IMPORTANT NOTE ABOUT PARALLELIZATION" above
## in this man page.

setAutoRealizationBackend()

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

grid <- defaultAutoGrid(sink, block.length=600)
for (bid in seq_along(grid)) {
    viewport <- grid[[bid]]
    block <- 101 * bid + runif(length(viewport))
    dim(block) <- dim(viewport)
    sink <- write_block(sink, viewport, block)
}

close(sink)
A <- as(sink, "DelayedArray")
A
}
\keyword{utilities}
