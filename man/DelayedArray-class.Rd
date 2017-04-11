\name{DelayedArray-class}
\docType{class}

\alias{class:DelayedArray}
\alias{DelayedArray-class}
\alias{DelayedArray}

\alias{class:DelayedMatrix}
\alias{DelayedMatrix-class}
\alias{DelayedMatrix}

\alias{coerce,DelayedArray,DelayedMatrix-method}

\alias{dim,DelayedArray-method}
\alias{dim<-,DelayedArray-method}
\alias{length,DelayedArray-method}
\alias{isEmpty,DelayedArray-method}
\alias{dimnames,DelayedArray-method}
\alias{dimnames<-,DelayedArray-method}
\alias{names,DelayedArray-method}
\alias{names<-,DelayedArray-method}

\alias{drop,DelayedArray-method}

\alias{[,DelayedArray-method}
\alias{t,DelayedArray-method}

\alias{as.array.DelayedArray}
\alias{as.array,DelayedArray-method}
\alias{as.vector.DelayedArray}
\alias{as.vector,DelayedArray-method}
\alias{as.matrix.DelayedArray}
\alias{as.matrix,DelayedArray-method}
\alias{as.data.frame.DelayedArray}
\alias{as.data.frame,DelayedArray-method}

\alias{type}
\alias{type,array-method}
\alias{type,DelayedArray-method}

\alias{[[,DelayedArray-method}

\alias{show,DelayedArray-method}

\alias{c,DelayedArray-method}
\alias{splitAsList,DelayedArray-method}
\alias{split.DelayedArray}
\alias{split,DelayedArray,ANY-method}

% Internal stuff
\alias{matrixClass}
\alias{matrixClass,DelayedArray-method}

\alias{subset_seed_as_array}
\alias{subset_seed_as_array,ANY-method}
\alias{subset_seed_as_array,array-method}
\alias{subset_seed_as_array,data.frame-method}
\alias{subset_seed_as_array,DataFrame-method}

\title{DelayedArray objects}

\description{
  Wrapping an array-like object (typically an on-disk object) in a
  DelayedArray object allows one to perform common array operations on it
  without loading the object in memory. In order to reduce memory usage and
  optimize performance, operations on the object are either delayed or
  executed using a block processing mechanism.
}

\usage{
DelayedArray(seed)  # constructor function
}

\arguments{
  \item{seed}{
    An array-like object.
  }
}

\section{In-memory versus on-disk realization}{
  To \emph{realize} a DelayedArray object (i.e. to trigger execution of the
  delayed operations carried by the object and return the result as an
  ordinary array), call \code{as.array} on it. However this realizes the
  full object at once \emph{in memory} which could require too much memory
  if the object is big. A big DelayedArray object is preferrably realized
  \emph{on disk} e.g. by calling \code{\link[HDF5Array]{writeHDF5Array}} on
  it (this function is defined in the \pkg{HDF5Array} package) or coercing it
  to an \link[HDF5Array]{HDF5Array} object with \code{as(x, "HDF5Array")}.
  Other on-disk backends can be supported. This uses a block-processing
  strategy so that the full object is not realized at once in memory. Instead
  the object is processed block by block i.e. the blocks are realized in
  memory and written to disk one at a time.
  See \code{?\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array} package
  for more information about this.
}

\section{Accessors}{
  DelayedArray objects support the same set of getters as ordinary arrays
  i.e. \code{dim()}, \code{length()}, and \code{dimnames()}.

  Only \code{dimnames()} is supported as a setter.
}

\section{Subsetting}{
  A DelayedArray object can be subsetted like an ordinary object but with
  the following differences:
  \itemize{
    \item The \code{drop} argument of the \code{[} operator is ignored i.e.
          subsetting a DelayedArray object always returns a DelayedArray
          object with the same number of dimensions. You need to call
          \code{drop()} on the subsetted object to actually drop its
          ineffective dimensions (i.e. the dimensions equal to 1).

    \item Linear subsetting (a.k.a. 1D-style subsetting, that is, subsetting
          with a single subscript \code{i}) is not supported.
  }

  Subsetting with \code{[[} is supported but only the linear form of it.

  DelayedArray objects don't support subassignment (\code{[<-} or \code{[[<-}).
}

\seealso{
  \itemize{
    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link[DelayedArray]{cbind} in this package (\pkg{DelayedArray})
          for binding DelayedArray objects along their rows or columns.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writting an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created 
          HDF5 datasets.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## WRAP AN ORDINARY ARRAY IN A DelayedArray OBJECT
## ---------------------------------------------------------------------
a <- array(runif(1500000), c(10000, 30, 5))
A <- DelayedArray(a)
A

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B             # still 3 dimensions (subsetting a DelayedArray object
              # never drops dimensions)
B <- drop(B)
B

cs <- colSums(b)
CS <- colSums(B)
stopifnot(identical(cs, CS))

## ---------------------------------------------------------------------
## WRAP A DataFrame OBJECT IN A DelayedArray OBJECT
## ---------------------------------------------------------------------

## Generate random coverage and score along an imaginary chromosome:
cov <- Rle(sample(20, 5000, replace=TRUE), sample(6, 5000, replace=TRUE))
score <- Rle(sample(100, nrun(cov), replace=TRUE), runLength(cov))

DF <- DataFrame(cov, score)
A2 <- DelayedArray(DF)
A2

t(A2)  # delayed transpose is very fast and very memory efficient because
       # the matrix data is not copied
colSums(A2)

## ---------------------------------------------------------------------
## A HDF5Array OBJECT IS A (PARTICULAR KIND OF) DelayedArray OBJECT
## ---------------------------------------------------------------------
library(HDF5Array)
A3 <- as(a, "HDF5Array")   # write 'a' to an HDF5 file
A3
is(A3, "DelayedArray")     # TRUE

B3 <- toto(A3)             # very fast! (operations are delayed)

B3                         # not a HDF5Array object because now it
                           # carries delayed operations
B3 <- drop(B3)

CS3 <- colSums(B3)
stopifnot(identical(cs, CS3))

## ---------------------------------------------------------------------
## PERFORM THE DELAYED OPERATIONS
## ---------------------------------------------------------------------
as(B3, "HDF5Array")        # "realize" 'B3' on disk

## If this is just an intermediate result, you can either keep going
## with B3 or replace it with its "realized" version:
B3 <- as(B3, "HDF5Array")  # no more delayed operations on new 'B3'

## For convenience, the backend used for realization can be set at the
## session level and realize() used instead of explicit coercion:
D <- cbind(B3, exp(B3))
setRealizeBackend("HDF5Array")
D <- realize(D)
}
\keyword{methods}
\keyword{classes}
