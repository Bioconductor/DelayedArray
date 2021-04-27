\name{AutoBlock-global-settings}

\alias{AutoBlock-global-settings}

\alias{getAutoBlockSize}
\alias{setAutoBlockSize}
\alias{get_type_size}
\alias{getAutoBlockLength}
\alias{getAutoBlockShape}
\alias{setAutoBlockShape}

\title{Control the geometry of automatic blocks}

\description{
  A family of utilities to control the automatic block size (or
  length) and shape.
}

\usage{
getAutoBlockSize()
setAutoBlockSize(size=1e8)

getAutoBlockLength(type)

getAutoBlockShape()
setAutoBlockShape(shape=c("hypercube",
                          "scale",
                          "first-dim-grows-first",
                          "last-dim-grows-first"))
}

\arguments{
  \item{size}{
    The \emph{auto block size} (automatic block size) in bytes. Note that,
    except when the type of the array data is \code{"character"} or
    \code{"list"}, the size of a block is its length multiplied by the
    size of an array element. For example, a block of 500 x 1000 x 500
    doubles has a length of 250 million elements and a size of 2 Gb (each
    double occupies 8 bytes of memory).

    The \emph{auto block size} is set to 100 Mb at package startup and can
    be reset anytime to this value by calling \code{setAutoBlockSize()}
    with no argument.
  }
  \item{type}{
    A string specifying the type of the array data.
  }
  \item{shape}{
    A string specifying the \emph{auto block shape} (automatic block shape).
    See \code{\link{makeCappedVolumeBox}} for a description of the
    supported shapes.

    The \emph{auto block shape} is set to \code{"hypercube"} at
    package startup and can be reset anytime to this value by calling
    \code{setAutoBlockShape()} with no argument.
  }
}

\details{
  \emph{block size} != \emph{block length}

  \emph{block length} = number of array elements in a block
  (i.e. \code{prod(dim(block))}).

  \emph{block size} = \emph{block length} * size of the individual elements
  in memory.

  For example, for an integer array, \emph{block size} (in bytes) is
  going to be 4 x \emph{block length}. For a numeric array \code{x}
  (i.e. \code{type(x) == "double"}), it's going to be 8 x \emph{block length}.

  In its current form, block processing in the \pkg{DelayedArray} package
  must decide the geometry of the blocks before starting the walk on the
  blocks. It does this based on several criteria. Two of them are:
  \itemize{
    \item The \emph{auto block size}: maximum size (in bytes) of a block
          once loaded in memory.
    \item The \code{type()} of the array (e.g. \code{integer}, \code{double},
          \code{complex}, etc...)
  }

  The \emph{auto block size} setting and \code{type(x)} control the maximum
  length of the blocks. Other criteria control their shape. So for example
  if you set the \emph{auto block size} to 8GB, this will cap the length of
  the blocks to 2e9 if your DelayedArray object \code{x} is of type
  \code{integer}, and to 1e9 if it's of type \code{double}.

  Note that this simple relationship between \emph{block size} and
  \emph{block length} assumes that blocks are loaded in memory as
  ordinary (a.k.a. dense) matrices or arrays. With sparse blocks,
  all bets are off. But the max block length is always taken to be
  the \emph{auto block size} divided by \code{get_type_size(type())}
  whether the blocks are going to be loaded as dense or sparse arrays.
  If they are going to be loaded as sparse arrays, their memory footprint
  is very likely to be smaller than if they were loaded as dense arrays
  so this is safe (although probably not optimal).

  It's important to keep in mind that the \emph{auto block size} setting
  is a simple way for the user to put a cap on the memory footprint of
  the blocks. Nothing more. In particular it doesn't control the maximum
  amount of memory used by the block processing algorithm. Other variables
  can impact dramatically memory usage like parallelization (where more than
  one block is loaded in memory at any given time), what the algorithm is
  doing with the blocks (e.g. something like \code{blockApply(x, identity)}
  will actually load the entire array data in memory), what delayed
  operations are on \code{x}, etc... It would be awesome to have a way to
  control the maximum amount of memory used by a block processing algorithm
  as a whole but we don't know how to do that.
}

\value{
  \code{getAutoBlockSize}: The current \emph{auto block size} in bytes
  as a single numeric value.

  \code{setAutoBlockSize}: The new \emph{auto block size} in bytes as an
  invisible single numeric value.

  \code{getAutoBlockLength}: The \emph{auto block length} as a single
  integer value.

  \code{getAutoBlockShape}: The current \emph{auto block shape} as a
  single string.

  \code{setAutoBlockShape}: The new \emph{auto block shape} as an invisible
  single string.
}

\seealso{
  \itemize{
    \item \code{\link{defaultAutoGrid}} and family to create automatic
          grids to use for block processing of array-like objects.

    \item \code{\link{blockApply}} and family for convenient block
          processing of an array-like object.

    \item The \code{\link{makeCappedVolumeBox}} utility to make
          \emph{capped volume boxes}.
  }
}

\examples{
getAutoBlockSize()

getAutoBlockLength("double")
getAutoBlockLength("integer")
getAutoBlockLength("logical")
getAutoBlockLength("raw")

m <- matrix(runif(600), ncol=12)
setAutoBlockSize(140)
getAutoBlockLength(type(m))
defaultAutoGrid(m)
lengths(defaultAutoGrid(m))
dims(defaultAutoGrid(m))

getAutoBlockShape()
setAutoBlockShape("scale")
defaultAutoGrid(m)
lengths(defaultAutoGrid(m))
dims(defaultAutoGrid(m))

## Reset the auto block size and shape to factory settings:
setAutoBlockSize()
setAutoBlockShape()
}
\keyword{utilities}
