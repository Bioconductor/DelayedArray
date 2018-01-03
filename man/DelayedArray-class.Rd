\name{DelayedArray-class}
\docType{class}

\alias{class:DelayedArray}
\alias{DelayedArray-class}

\alias{class:DelayedMatrix}
\alias{DelayedMatrix-class}
\alias{DelayedMatrix}

\alias{coerce,DelayedArray,DelayedMatrix-method}

\alias{DelayedArray}
\alias{DelayedArray,ANY-method}
\alias{DelayedArray,DelayedArray-method}

\alias{updateObject,DelayedArray-method}

\alias{seed}
\alias{seed,DelayedArray-method}
\alias{seed<-}
\alias{seed<-,DelayedArray-method}

\alias{path}
\alias{path,DelayedArray-method}
\alias{path<-,DelayedArray-method}

\alias{type}

\alias{dim,DelayedArray-method}
\alias{isEmpty,DelayedArray-method}

\alias{coerce,DelayedArray,SeedDimCombiner-method}
\alias{aperm}
\alias{aperm.DelayedArray}
\alias{aperm,DelayedArray-method}

\alias{dim<-,DelayedArray-method}
\alias{dimnames,DelayedArray-method}
\alias{dimnames<-,DelayedArray-method}
\alias{names,DelayedArray-method}
\alias{names<-,DelayedArray-method}

\alias{drop,DelayedArray-method}

\alias{[,DelayedArray-method}
\alias{[<-,DelayedArray-method}

\alias{as.array.DelayedArray}
\alias{as.array,DelayedArray-method}
\alias{as.matrix.DelayedArray}
\alias{as.matrix,DelayedArray-method}
\alias{as.data.frame.DelayedArray}
\alias{as.data.frame,DelayedArray-method}
\alias{as.vector.DelayedArray}
\alias{as.vector,DelayedArray-method}
\alias{as.logical.DelayedArray}
\alias{as.logical,DelayedArray-method}
\alias{as.integer.DelayedArray}
\alias{as.integer,DelayedArray-method}
\alias{as.numeric.DelayedArray}
\alias{as.numeric,DelayedArray-method}
\alias{as.complex.DelayedArray}
\alias{as.complex,DelayedArray-method}
\alias{as.character.DelayedArray}
\alias{as.character,DelayedArray-method}
\alias{as.raw.DelayedArray}
\alias{as.raw,DelayedArray-method}

\alias{coerce,DelayedMatrix,dgCMatrix-method}
\alias{coerce,DelayedMatrix,sparseMatrix-method}

\alias{[[,DelayedArray-method}

\alias{show,DelayedArray-method}

\alias{c,DelayedArray-method}
\alias{splitAsList,DelayedArray-method}
\alias{split.DelayedArray}
\alias{split,DelayedArray,ANY-method}

\alias{rbind}
\alias{rbind,DelayedMatrix-method}
\alias{rbind,DelayedArray-method}
\alias{arbind,DelayedArray-method}

\alias{cbind}
\alias{cbind,DelayedMatrix-method}
\alias{cbind,DelayedArray-method}
\alias{acbind,DelayedArray-method}

% Internal stuff
\alias{matrixClass}
\alias{matrixClass,DelayedArray-method}

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
seed(x)             # seed getter
path(object, ...)   # path getter
type(x)
}

\arguments{
  \item{seed}{
    An array-like object.
  }
  \item{x, object}{
    A DelayedArray object.
  }
  \item{...}{
    Additional arguments passed to methods.
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
  In addition, they support \code{seed()}, \code{path()}, and \code{type()}.
  \code{type()} is the DelayedArray equivalent of \code{typeof()} (or
  \code{storage.mode()}) for ordinary arrays. Note that, for convenience
  and consistency, \code{type()} also works on ordinary arrays.

  \code{dimnames()}, \code{seed()}, and \code{path()} also work as setters.
}

\section{Subsetting}{
  A DelayedArray object can be subsetted with \code{[} like an ordinary array
  but with the following differences:
  \itemize{
    \item \emph{Multi-dimensional single bracket subsetting} (i.e. subsetting
          of the form \code{x[i_1, i_2, ..., i_n]} with one (possibly missing)
          subscript per dimension) returns a DelayedArray object where the
          subsetting is actually delayed. So it's a very light operation.

    \item The \code{drop} argument of the \code{[} operator is ignored i.e.
          subsetting a DelayedArray object always returns a DelayedArray
          object with the same number of dimensions as the original object.
          You need to call \code{drop()} on the subsetted object to actually
          drop its ineffective dimensions (i.e. the dimensions equal to 1).
          \code{drop()} is also a delayed operation so is very light.

    \item \emph{Linear single bracket subsetting} (a.k.a. 1D-style subsetting,
          that is, subsetting of the form \code{x[i]}) only works if subscript
          \code{i} is a numeric vector at the moment. Furthermore, \code{i}
          cannot contain NAs and all the indices in it must be >= 1 and <=
          \code{length(x)} for now. It returns an atomic vector of the same
          length as \code{i}. This is NOT a delayed operation.
  }

  Subsetting with \code{[[} is supported but only the \emph{linear} form
  of it at the moment i.e. the \code{x[[i]]} form where \code{i} is a
  \emph{single} numeric value >= 1 and <= \code{length(x)}. It is equivalent
  to \code{x[i]}.

  DelayedArray objects support only 2 forms of subassignment at the moment:
  \code{x[i] <- value} and \code{x[] <- value}. The former is supported only
  when the subscript \code{i} is a logical DelayedArray object with the same
  dimensions as \code{x} and when \code{value} is a \emph{scalar} (i.e. an
  atomic vector of length 1). The latter is supported only when \code{value}
  is an atomic vector and \code{length(value)} is a divisor of \code{nrow(x)}.
  Both are delayed operations so are very light.

  Single value replacement (\code{x[[...]] <- value}) is not supported.
}

\section{Binding}{
  Binding DelayedArray objects along the rows (or columns) is supported
  via the \code{rbind} and \code{arbind} (or \code{cbind} and \code{acbind})
  methods for DelayedArray objects. All these operations are delayed.
}

\seealso{
  \itemize{
    \item \code{\link{type}} to get the type of the elements of an array-like
          object.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \code{\link[base]{cbind}} in the \pkg{base} package for
          rbind/cbind'ing ordinary arrays.

    \item \code{\link{acbind}} in this package (\pkg{DelayedArray}) for
          arbind/acbind'ing ordinary arrays.

    \item \link{RleArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. WRAP AN ORDINARY ARRAY IN A DelayedArray OBJECT
## ---------------------------------------------------------------------
a <- array(runif(1500000), dim=c(10000, 30, 5))
A <- DelayedArray(a)
A
## The seed of A is treated as a "read-only" object so won't change when
## we start operating on A:
stopifnot(identical(a, seed(A)))
type(A)

## Multi-dimensional single bracket subsetting:
m <- a[11:20 , 5, ]        # a matrix
A[11:20 , 5, ]             # not a DelayedMatrix (still 3 dimensions)
M <- drop(A[11:20 , 5, ])  # a DelayedMatrix object
stopifnot(identical(m, as.array(M)))

## Linear single bracket subsetting:
A[11:20]
A[which(A <= 1e-5)]

## Subassignment:
A[A < 0.2] <- NA
a[a < 0.2] <- NA
stopifnot(identical(a, as.array(A)))

## Other operations:
toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B  # still 3 dimensions (subsetting a DelayedArray object never drops
   # dimensions)
B <- drop(B)
B

cs <- colSums(b)
CS <- colSums(B)
stopifnot(identical(cs, CS))

## ---------------------------------------------------------------------
## B. WRAP A DataFrame OBJECT IN A DelayedArray OBJECT
## ---------------------------------------------------------------------

## Generate random coverage and score along an imaginary chromosome:
cov <- Rle(sample(20, 5000, replace=TRUE), sample(6, 5000, replace=TRUE))
score <- Rle(sample(100, nrun(cov), replace=TRUE), runLength(cov))

DF <- DataFrame(cov, score)
A2 <- DelayedArray(DF)
A2
seed(A2)  # 'DF'

## Coercion of a DelayedMatrix object to DataFrame produces a DataFrame
## object with Rle columns:
as(A2, "DataFrame")
stopifnot(identical(DF, as(A2, "DataFrame")))

t(A2)  # transposition is delayed so is very fast and memory efficient
stopifnot(identical(A2, t(t(A2))))
colSums(A2)

## ---------------------------------------------------------------------
## C. A HDF5Array OBJECT IS A (PARTICULAR KIND OF) DelayedArray OBJECT
## ---------------------------------------------------------------------
library(HDF5Array)
A3 <- as(a, "HDF5Array")   # write 'a' to an HDF5 file
A3
is(A3, "DelayedArray")     # TRUE
seed(A3)                   # an HDF5ArraySeed object

B3 <- toto(A3)             # very fast! (operations are delayed)

B3                         # not a HDF5Array object anymore because now
                           # it carries delayed operations
B3 <- drop(B3)

CS3 <- colSums(B3)
stopifnot(identical(cs, CS3))

## ---------------------------------------------------------------------
## D. PERFORM THE DELAYED OPERATIONS
## ---------------------------------------------------------------------
as(B3, "HDF5Array")        # "realize" 'B3' on disk

## If this is just an intermediate result, you can either keep going
## with B3 or replace it with its "realized" version:
B3 <- as(B3, "HDF5Array")  # no more delayed operations on new 'B3'
seed(B3)
path(B3)

## For convenience, realize() can be used instead of explicit coercion.
## The current "realization backend" controls where realization
## happens e.g. in memory if set to NULL or in an HDF5 file if set
## to "HDF5Array":
D <- cbind(B3, exp(B3))
D
setRealizationBackend("HDF5Array")
D <- realize(D)
D
## See '?realize' for more information about "realization backends".

## ---------------------------------------------------------------------
## E. BIND DelayedArray OBJECTS
## ---------------------------------------------------------------------

## rbind/cbind

library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M12 <- rbind(M1, t(M2))
M12
colMeans(M12)

## arbind/acbind

example(acbind)  # to create arrays a1, a2, a3

A1 <- DelayedArray(a1)
A2 <- DelayedArray(a2)
A3 <- DelayedArray(a3)
A <- arbind(A1, A2, A3)
A

## Sanity check:
stopifnot(identical(arbind(a1, a2, a3), as.array(A)))

## ---------------------------------------------------------------------
## F. MODIFY THE PATH OF A DelayedArray OBJECT
## ---------------------------------------------------------------------
## This can be useful if the file containing the array data is on a
## shared partition but the exact path to the partition depends on the
## machine from which the data is being accessed.
## For example:

\dontrun{
A <- HDF5Array("/path/to/lab_data/my_precious_data.h5")
path(A)

## Operate on A...
## Now A carries delayed operations.
## Make sure path(A) still works:
path(A)

## Save A:
save(A, file="A.rda")

## A.rda should be small (it doesn't contain the array data).
## Send it to a co-worker that has access to my_precious_data.h5.

## Co-worker loads it:
load("A.rda")
path(A)

## A is broken because path(A) is incorrect for co-worker:
A  # error!

## Co-worker fixes the path (in this case this is better done using the
## dirname() setter rather than the path() setter):
dirname(A) <- "E:/other/path/to/lab_data"

## A "works" again:
A
}

## ---------------------------------------------------------------------
## G. WRAP A SPARSE MATRIX IN A DelayedArray OBJECT
## ---------------------------------------------------------------------
\dontrun{
library(Matrix)
M <- 75000L
N <- 1800L
p <- sparseMatrix(sample(M, 9000000, replace=TRUE),
                  sample(N, 9000000, replace=TRUE),
                  x=runif(9000000), dims=c(M, N))
P <- DelayedArray(p)
P
p2 <- as(P, "sparseMatrix")
stopifnot(identical(p, p2))

## The following is based on the following post by Murat Tasan on the
## R-help mailing list:
##   https://stat.ethz.ch/pipermail/r-help/2017-May/446702.html

## As pointed out by Murat, the straight-forward row normalization
## directly on sparse matrix 'p' would consume too much memory:
row_normalized_p <- p / rowSums(p^2)  # consumes too much memory
## because the rowSums() result is being recycled (appropriately) into a
## *dense* matrix with dimensions equal to dim(p).

## Murat came up with the following solution that is very fast and memory
## efficient:
row_normalized_p1 <- Diagonal(x=1/sqrt(Matrix::rowSums(p^2))) %*% p

## With a DelayedArray object, the straight-forward approach uses a
## block processing strategy behind the scene so it doesn't consume
## too much memory.

## First, let's see the block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
## and set block size to a bigger value than the default:
getOption("DelayedArray.block.size")
options(DelayedArray.block.size=80e6)

row_normalized_P <- P / sqrt(DelayedArray::rowSums(P^2))

## Increasing the block size increases the speed but also memory usage:
options(DelayedArray.block.size=200e6)
row_normalized_P2 <- P / sqrt(DelayedArray::rowSums(P^2))
stopifnot(all.equal(row_normalized_P, row_normalized_P2))

## Back to sparse representation:
DelayedArray:::set_verbose_block_processing(FALSE)
row_normalized_p2 <- as(row_normalized_P, "sparseMatrix")
stopifnot(all.equal(row_normalized_p1, row_normalized_p2))

options(DelayedArray.block.size=10e6)
}
}
\keyword{methods}
\keyword{classes}
