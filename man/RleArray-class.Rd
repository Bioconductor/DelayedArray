\name{RleArray-class}
\docType{class}

\alias{class:RleArraySeed}
\alias{RleArraySeed-class}
\alias{class:SolidRleArraySeed}
\alias{SolidRleArraySeed-class}

\alias{class:RleRealizationSink}
\alias{RleRealizationSink-class}
\alias{type,RleRealizationSink-method}

\alias{class:ChunkedRleArraySeed}
\alias{ChunkedRleArraySeed-class}

\alias{dim,RleArraySeed-method}
\alias{dimnames,RleArraySeed-method}

\alias{coerce,SolidRleArraySeed,Rle-method}
\alias{coerce,RleRealizationSink,RleList-method}
\alias{coerce,RleRealizationSink,Rle-method}

\alias{extract_array,SolidRleArraySeed-method}
\alias{extract_array,ChunkedRleArraySeed-method}

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

\alias{DelayedArray,RleArraySeed-method}
\alias{RleArray}
\alias{coerce,RleArray,Rle-method}

\alias{write_block,RleRealizationSink-method}
\alias{coerce,RleRealizationSink,RleArray-method}
\alias{coerce,RleRealizationSink,DelayedArray-method}
\alias{coerce,ANY,RleArray-method}
\alias{coerce,DelayedArray,RleArray-method}
\alias{coerce,DelayedMatrix,RleMatrix-method}

\alias{coerce,RleList,RleArray-method}
\alias{coerce,RleMatrix,RleList-method}

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
RleArray(data, dim, dimnames=NULL, chunksize=NULL)  # constructor function
}

\arguments{
  \item{data}{
    An \link[S4Vectors]{Rle} object, or an ordinary list of Rle objects,
    or an \link[IRanges]{RleList} object, or a \link[S4Vectors]{DataFrame}
    object where all the columns are Rle objects. More generally speaking,
    \code{data} can be any list-like object where all the list elements
    are Rle objects.
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
    \item \link[S4Vectors]{Rle} objects in the \pkg{S4Vectors} package and
          \link[IRanges]{RleList} objects in the \pkg{IRanges} package.

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
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

rle <- Rle(sample(6L, 500000, replace=TRUE), 8)
a <- array(rle, dim=c(50, 20, 4000))  # array() expands the Rle object
                                      # internally with as.vector()

A <- RleArray(rle, dim=c(50, 20, 4000))  # Rle object is NOT expanded
A

object.size(a)
object.size(A)

stopifnot(identical(a, as.array(A)))

as(A, "Rle")  # deconstruction

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
m1 <- toto(a)
head(m1)

M1 <- toto(A)  # very fast! (operations are delayed)
M1

stopifnot(identical(m1, as.array(M1)))

cs <- colSums(m1)
CS <- colSums(M1)
stopifnot(identical(cs, CS))

## Coercing a DelayedMatrix object to DataFrame produces a DataFrame
## object with Rle columns:
as(M1, "DataFrame")

## Coercing an RleList object to RleMatrix only works if all the list
## elements in the former have the same length. It returns an RleMatrix
## object with 1 column per list element in the RleList object. The
## names on the RleList object are propagated as the colnames of the
## returned object.
rle_list <- RleList(A=Rle(sample(3L, 100, replace=TRUE)),
                    B=Rle(sample(3L, 100, replace=TRUE)))
M2 <- as(rle_list, "RleMatrix")
M2
stopifnot(identical(as(M2, "RleList"), rle_list))

## ---------------------------------------------------------------------
## CONSTRUCTING A LARGE RleArray OBJECT
## ---------------------------------------------------------------------

## The RleArray() constructor does not accept a long Rle object at the
## moment:
\dontrun{
  RleArray(Rle(5, 3e9), dim=c(3, 1e9))  # error!
}

## The workaround is to supply a list of Rle objects instead:
data <- lapply(1:500, function(j) Rle(runif(99), 1e6 + 99:1))
dim <- c(6750, 73337, 100)
A <- RleArray(data, dim)
A

## If we first turn 'data' into an RleList object, then we can use
## this object to construct an RleMatrix where each column corresponds
## to a list element in the RleList object:
rle_list <- as(data, "RleList")
M3 <- as(rle_list, "RleMatrix")
M3
stopifnot(identical(as(M3, "RleList"), rle_list))
}
\keyword{classes}
\keyword{methods}
