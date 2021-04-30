\name{RleArray-class}
\docType{class}

\alias{class:RleArray}
\alias{RleArray-class}
\alias{RleArray}

\alias{DelayedArray,RleArraySeed-method}
\alias{coerce,RleArray,Rle-method}

\alias{class:RleMatrix}
\alias{RleMatrix-class}
\alias{RleMatrix}

\alias{matrixClass,RleArray-method}
\alias{coerce,RleArray,RleMatrix-method}
\alias{coerce,RleMatrix,RleArray-method}
\alias{coerce,ANY,RleMatrix-method}

\alias{coerce,RleList,RleArray-method}
\alias{coerce,RleMatrix,RleList-method}

\alias{coerce,DataFrame,RleArray-method}
\alias{coerce,RleMatrix,DataFrame-method}
\alias{coerce,DelayedMatrix,DataFrame-method}

\alias{write_block,RleRealizationSink-method}
\alias{coerce,RleRealizationSink,RleArray-method}
\alias{coerce,RleRealizationSink,DelayedArray-method}
\alias{coerce,ANY,RleArray-method}
\alias{coerce,DelayedArray,RleArray-method}
\alias{coerce,DelayedMatrix,RleMatrix-method}

\title{RleArray objects}

\description{
  The RleArray class is a \link[DelayedArray]{DelayedArray} subclass
  for representing an in-memory Run Length Encoded array-like dataset.

  All the operations available for \link[DelayedArray]{DelayedArray}
  objects work on RleArray objects.
}

\usage{
## Constructor function:
RleArray(data, dim, dimnames, chunksize=NULL)
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
    The \emph{dimnames} of the object to be created. Must be \code{NULL} or
    a list of length the number of dimensions. Each list element must be
    either \code{NULL} or a character vector along the corresponding dimension.
  }
  \item{chunksize}{
    Experimental. Don't use!
  }
}

\value{
  An RleArray (or RleMatrix) object. (Note that RleMatrix extends RleArray.)
}

\seealso{
  \itemize{
    \item \link[S4Vectors]{Rle} and \link[S4Vectors]{DataFrame} objects
          in the \pkg{S4Vectors} package and \link[IRanges]{RleList} objects
          in the \pkg{IRanges} package.

    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.

    \item \link{ConstantArray} objects for mimicking an array containing
          a constant value, without actually creating said array in memory.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item The \link{RleArraySeed} helper class.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BASIC EXAMPLE
## ---------------------------------------------------------------------

data <- Rle(sample(6L, 500000, replace=TRUE), 8)
a <- array(data, dim=c(50, 20, 4000))  # array() expands the Rle object
                                       # internally with as.vector()

A <- RleArray(data, dim=c(50, 20, 4000))  # Rle object is NOT expanded
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

## ---------------------------------------------------------------------
## B. MAKING AN RleArray OBJECT FROM A LIST-LIKE OBJECT OF Rle OBJECTS
## ---------------------------------------------------------------------

## From a DataFrame object:
DF <- DataFrame(A=Rle(sample(3L, 100, replace=TRUE)),
                B=Rle(sample(3L, 100, replace=TRUE)),
                C=Rle(sample(3L, 100, replace=TRUE) - 0.5),
                row.names=sprintf("ID\%03d", 1:100))

M2 <- RleArray(DF)
M2

A3 <- RleArray(DF, dim=c(25, 6, 2))
A3

M4 <- RleArray(DF, dim=c(25, 12), dimnames=list(LETTERS[1:25], NULL))
M4

## From an ordinary list:
## If all the supplied Rle objects have the same length and if the 'dim'
## argument is not specified, then the RleArray() constructor returns an
## RleMatrix object with 1 column per Rle object. If the 'dimnames'
## argument is not specified, then the names on the list are propagated
## as the colnames of the returned object.
data <- as.list(DF)
M2b <- RleArray(data)
A3b <- RleArray(data, dim=c(25, 6, 2))
M4b <- RleArray(data, dim=c(25, 12), dimnames=list(LETTERS[1:25], NULL))

data2 <- list(Rle(sample(3L, 9, replace=TRUE)) * 11L,
              Rle(sample(3L, 15, replace=TRUE)))
\dontrun{
  RleArray(data2)  # error! (cannot infer the dim)
}
RleArray(data2, dim=c(4, 6))

## From an RleList object:
data <- RleList(data)
M2c <- RleArray(data)
A3c <- RleArray(data, dim=c(25, 6, 2))
M4c <- RleArray(data, dim=c(25, 12), dimnames=list(LETTERS[1:25], NULL))

data2 <- RleList(data2)
\dontrun{
  RleArray(data2)  # error! (cannot infer the dim)
}
RleArray(data2, dim=4:2)

## Sanity checks:
data0 <- as.vector(unlist(DF, use.names=FALSE))
m2 <- matrix(data0, ncol=3, dimnames=dimnames(M2))
stopifnot(identical(m2, as.matrix(M2)))
rownames(m2) <- NULL
stopifnot(identical(m2, as.matrix(M2b)))
stopifnot(identical(m2, as.matrix(M2c)))
a3 <- array(data0, dim=c(25, 6, 2))
stopifnot(identical(a3, as.array(A3)))
stopifnot(identical(a3, as.array(A3b)))
stopifnot(identical(a3, as.array(A3c)))
m4 <- matrix(data0, ncol=12, dimnames=dimnames(M4))
stopifnot(identical(m4, as.matrix(M4)))
stopifnot(identical(m4, as.matrix(M4b)))
stopifnot(identical(m4, as.matrix(M4c)))

## ---------------------------------------------------------------------
## C. COERCING FROM RleList OR DataFrame TO RleMatrix
## ---------------------------------------------------------------------

## Coercing an RleList object to RleMatrix only works if all the list
## elements in the former have the same length.
x <- RleList(A=Rle(sample(3L, 20, replace=TRUE)),
             B=Rle(sample(3L, 20, replace=TRUE)))
M <- as(x, "RleMatrix")
stopifnot(identical(x, as(M, "RleList")))

x <- DataFrame(A=x[[1]], B=x[[2]], row.names=letters[1:20])
M <- as(x, "RleMatrix")
stopifnot(identical(x, as(M, "DataFrame")))

## ---------------------------------------------------------------------
## D. CONSTRUCTING A LARGE RleArray OBJECT
## ---------------------------------------------------------------------

## The RleArray() constructor does not accept a "long" Rle object (i.e.
## an object of length > .Machine$integer.max) at the moment:
\dontrun{
  RleArray(Rle(5, 3e9), dim=c(3, 1e9))  # error!
}

## The workaround is to supply a list of Rle objects instead:

toy_Rle <- function() {
  run_lens <- c(sample(4), sample(rep(c(1:19, 40) * 3, 6e4)), sample(4))
  run_vals <- sample(700, length(run_lens), replace=TRUE) / 5
  Rle(run_vals, run_lens)
}
rle_list <- lapply(1:80, function(j) toy_Rle())  # takes about 20 sec.

## Cumulative length of all the Rle objects is > .Machine$integer.max:
sum(lengths(rle_list))  # 3.31e+09

## Feed 'rle_list' to the RleArray() constructor:
dim <- c(14395, 320, 719)
A <- RleArray(rle_list, dim)
A

## Because all the Rle objects in 'rle_list' have the same length, we
## can call RleArray() on it without specifying the 'dim' argument. This
## returns an RleMatrix object where each column corresponds to an Rle
## object in 'rle_list':
M <- RleArray(rle_list)
M
stopifnot(identical(as(rle_list, "RleList"), as(M, "RleList")))

## ---------------------------------------------------------------------
## E. CHANGING THE TYPE OF AN RleArray OBJECT FROM "double" TO "integer"
## ---------------------------------------------------------------------

## An RleArray object is an in-memory object so it can be useful to
## reduce its memory footprint. For an object of type "double" this can
## be done by changing its type to "integer" (integers are half the size
## of doubles in memory). Of course this only makes sense if this results
## in a loss of precision that is acceptable.
## On an ordinary array (or matrix) 'a', this is simply a matter of
## doing 'storage.mode(a) <- "integer"'. However, with a DelayedArray
## object, things are a little bit different. Let's do this on a subset
## of the RleMatrix object 'M' created in the previous section.

M1 <- as(M[1:6e5, ], "RleMatrix")
rm(M)

## First of all, it's important to be aware that object.size() (from
## package utils) is NOT reliable on RleArray objects! This is because
## the data in an RleArray object is stored in an environment and
## object.size() stubbornly refuses to take the content of an environment
## into account when computing its size:
object.size(list2env(list(aa=1:10)))   # 56 bytes
object.size(list2env(list(aa=1:1e6)))  # always 56 bytes!

## So we'll use obj_size() instead (from package lobstr):
library(lobstr)
obj_size(list2env(list(aa=1:10)))   # 264 B
obj_size(list2env(list(aa=1:1e6)))  # 4 MB
obj_size(list2env(list(aa=as.double(1:1e6))))  # 8 MB

obj_size(M1)  # 16.7 MB

type(M1) <- "integer"  # Delayed!
M1                     # Note the class: it's no longer RleMatrix!
                       # (That's because the object now carries delayed
                       # operations.)

## Because changing the type is a delayed operation, the memory footprint
## of the object has not changed yet (remember that the original data in
## a DelayedArray object is stored in its "seed" and its seed is never
## modified **in-place**, that is, no operation on the object will ever
## modify its seed):
obj_size(M1)  # Still the same (well, a very tiny more, because the
              # object is now carrying one more delayed operation,
              # the `type<-` operation)

## To effectively reduce the memory footprint of the object, a new object
## needs to be created. This is achieved simply by **realizing** M1 as a
## (new) RleArray object. Note that this realization will use block
## processing:

DelayedArray:::set_verbose_block_processing(TRUE)  # See block processing
                                                   # in action.
getAutoBlockSize()      # Automatic block size (100 Mb by default).
setAutoBlockSize(20e6)  # Set automatic block size to 20 Mb.

M2 <- as(M1, "RleArray")
DelayedArray:::set_verbose_block_processing(FALSE)
setAutoBlockSize()      # Reset automatic block size to factory settings.


M2
obj_size(M2)  # 6.91 MB (Less than half the original size! This is
              # because RleArray objects use some internal tricks to
              # reduce memory footprint even more when the data in
              # their seed is of type "integer".)

## Finally note that the 2-step approach described here (i.e.
## type(A) <- "integer" followed by realization) is generic and works
## on any kind of DelayedArray object or derivative. In particular,
## after doing 'type(A) <- "integer"', 'A' can be realized as anything
## as long as the realization backend is supported (e.g. could be
## 'as(A, "HDF5Array")' or 'as(A, "TENxMatrix")') and realization will
## always use block processing so the array data will never be fully
## loaded in memory.
}
\keyword{classes}
\keyword{methods}
