\name{RleArray-class}
\docType{class}

\alias{class:RleArraySeed}
\alias{RleArraySeed-class}
\alias{class:SolidRleArraySeed}
\alias{SolidRleArraySeed-class}
\alias{class:RleRealizationSink}
\alias{RleRealizationSink-class}
\alias{class:ChunkedRleArraySeed}
\alias{ChunkedRleArraySeed-class}

\alias{dim,RleArraySeed-method}
\alias{dimnames,RleArraySeed-method}

\alias{coerce,SolidRleArraySeed,Rle-method}
\alias{coerce,RleRealizationSink,Rle-method}

\alias{subset_seed_as_array,SolidRleArraySeed-method}
\alias{subset_seed_as_array,ChunkedRleArraySeed-method}

\alias{coerce,RleRealizationSink,ChunkedRleArraySeed-method}
\alias{coerce,ChunkedRleArraySeed,SolidRleArraySeed-method}

\alias{class:RleArray}
\alias{RleArray-class}

\alias{class:RleMatrix}
\alias{RleMatrix-class}
\alias{RleMatrix}

\alias{coerce,RleArray,RleMatrix-method}
\alias{matrixClass,RleArray-method}
\alias{coerce,ANY,RleMatrix-method}
\alias{coerce,RleArray,RleList-method}

\alias{DelayedArray,RleArraySeed-method}
\alias{RleArray}

\alias{write_block_to_sink,RleRealizationSink-method}
\alias{coerce,RleRealizationSink,RleArray-method}
\alias{coerce,RleRealizationSink,DelayedArray-method}
\alias{coerce,ANY,RleArray-method}
\alias{coerce,DelayedArray,RleArray-method}
\alias{coerce,DelayedMatrix,RleMatrix-method}

\alias{coerce,DataFrame,RleMatrix-method}
\alias{coerce,DataFrame,RleArray-method}
\alias{coerce,RleMatrix,DataFrame-method}
\alias{coerce,DelayedMatrix,DataFrame-method}

\title{RleArray objects}

\description{
  The RleArray class is an array-like container where the values are stored
  in a run-length encoding format. RleArray objects support delayed operations
  and block processing.
}

\usage{
RleArray(rle, dim, dimnames=NULL, chunksize=NULL)  # constructor function
}

\arguments{
  \item{rle}{
    An \link[S4Vectors]{Rle} object.
  }
  \item{dim}{
    The dimensions of the object to be created, that is, an integer vector
    of length one or more giving the maximal indices in each dimension.
  }
  \item{dimnames}{
    Either \code{NULL} or the names for the dimensions. This must a list of
    length the number of dimensions. Each list element must be either
    \code{NULL} or a character vector along the corresponding dimension.
  }
  \item{chunksize}{
    Experimental. Don't use!
  }
}

\details{
  RleArray extends \link{DelayedArray}. All the operations available on
  \link{DelayedArray} objects work on RleArray objects.
}

\seealso{
  \itemize{
    \item \link[S4Vectors]{Rle} objects in the \pkg{S4Vectors} package.

    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
rle <- Rle(sample(6L, 500000, replace=TRUE), 8)
a <- array(rle, dim=c(50, 20, 4000))  # array() expands the Rle object
                                      # internally with as.vector()

A <- RleArray(rle, dim=c(50, 20, 4000))  # Rle object is NOT expanded
A

object.size(a)
object.size(A)

stopifnot(identical(a, as.array(A)))

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B  # still 3 dimensions (subsetting a DelayedArray object never drops
   # dimensions)
B <- drop(B)
B

stopifnot(identical(b, as.array(B)))

cs <- colSums(b)
CS <- colSums(B)
stopifnot(identical(cs, CS))

## Coercion of a DelayedMatrix object to DataFrame produces a DataFrame
## object with Rle columns:
as(B, "DataFrame")

## Coercion of an RleList to an RleArray. Only works if each element of the
## RleList has the same length. Column names are taken from the names of the
## list elements
rle_list <- RleList(A = Rle(sample(3, 100, replace = TRUE)),
                    B = Rle(sample(3, 100, replace = TRUE)))
rle_array <- as(rle_list, "RleArray")
}
\keyword{methods}
\keyword{classes}
