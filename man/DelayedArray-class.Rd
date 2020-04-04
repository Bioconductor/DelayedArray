\name{DelayedArray-class}
\docType{class}

\alias{class:DelayedArray}
\alias{DelayedArray-class}

\alias{class:DelayedMatrix}
\alias{DelayedMatrix-class}
\alias{DelayedMatrix}

\alias{coerce,DelayedArray,DelayedMatrix-method}
\alias{coerce,DelayedMatrix,DelayedArray-method}

\alias{dim,DelayedArray-method}
\alias{dimnames,DelayedArray-method}
\alias{extract_array,DelayedArray-method}

\alias{new_DelayedArray}
\alias{DelayedArray}
\alias{DelayedArray,ANY-method}
\alias{DelayedArray,DelayedArray-method}
\alias{DelayedArray,DelayedOp-method}

\alias{class:DelayedArray1}
\alias{DelayedArray1-class}
\alias{DelayedArray1}
\alias{updateObject,DelayedArray-method}

\alias{type}

\alias{aperm}
\alias{aperm.DelayedArray}
\alias{aperm,DelayedArray-method}

\alias{dim<-,DelayedArray-method}
\alias{dimnames<-,DelayedArray-method}
\alias{names,DelayedArray-method}
\alias{names<-,DelayedArray-method}

\alias{drop,DelayedArray-method}
\alias{[,DelayedArray-method}
\alias{[<-,DelayedArray-method}

\alias{coerce,DelayedArray,SparseArraySeed-method}
\alias{coerce,DelayedMatrix,dgCMatrix-method}
\alias{coerce,DelayedMatrix,sparseMatrix-method}

\alias{[[,DelayedArray-method}

\alias{show,DelayedArray-method}

\alias{c,DelayedArray-method}
\alias{splitAsList,DelayedArray-method}
\alias{split.DelayedArray}
\alias{split,DelayedArray,ANY-method}

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
type(x)
}

\arguments{
  \item{seed}{
    An array-like object.
  }
  \item{x}{
    Typically a DelayedArray object. More generally \code{type()} is
    expected to work on any array-like object (that is, any object for
    which \code{dim(x)} is not NULL), or any ordinary vector (i.e.
    atomic or non-atomic).
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
  Other on-disk backends can be supported. This uses a block processing
  strategy so that the full object is not realized at once in memory. Instead
  the object is processed block by block i.e. the blocks are realized in
  memory and written to disk one at a time.
  See \code{?\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array} package
  for more information about this.
}

\section{Accessors}{
  DelayedArray objects support the same set of getters as ordinary arrays
  i.e. \code{dim()}, \code{length()}, and \code{dimnames()}.
  In addition, they support \code{type()}, \code{\link{nseed}()},
  \code{\link{seed}()}, and \code{\link{path}()}.

  \code{type()} is the DelayedArray equivalent of \code{typeof()} (or
  \code{storage.mode()}) for ordinary arrays and vectors. Note that, for
  convenience and consistency, \code{type()} also supports ordinary arrays
  and vectors. It should also support any array-like object, that is, any
  object \code{x} for which \code{dim(x)} is not NULL.

  \code{dimnames()}, \code{\link{seed}()}, and \code{\link{path}()} also
  work as setters.
}

\section{Subsetting}{
  A DelayedArray object can be subsetted with \code{[} like an ordinary array,
  but with the following differences:
  \itemize{
    \item \emph{Multi-dimensional single bracket subsetting} (i.e. subsetting
          of the form \code{x[i_1, i_2, ..., i_n]} with one (possibly missing)
          subscript per dimension) returns a DelayedArray object where the
          subsetting is actually delayed. So it's a very light operation.
          One notable exception to this is when \code{drop=TRUE} and the
          result has only one dimension, in which case it is returned as an
          ordinary vector (atomic or list).
          Note that NAs in the subscripts are not supported.

    \item \emph{Linear single bracket subsetting} (a.k.a. 1D-style subsetting,
          that is, subsetting of the form \code{x[i]}) only works if the
          subscript \code{i} is a numeric vector at the moment. Furthermore,
          \code{i} cannot contain NAs and all the indices in it must be >= 1
          and <= \code{length(x)} for now. It returns an atomic vector of the
          same length as \code{i}. This is NOT a delayed operation (block
          processing is triggered).
  }

  Subsetting with \code{[[} is supported but only the \emph{linear} form
  of it at the moment i.e. the \code{x[[i]]} form where \code{i} is a
  \emph{single} numeric value >= 1 and <= \code{length(x)}. It is equivalent
  to \code{x[i][[1]]}.

  Subassignment to a DelayedArray object with \code{[<-} is also supported
  like with an ordinary array, but with the following restrictions:
  \itemize{
    \item \emph{Multi-dimensional subassignment} (i.e. subassignment of the
          form \code{x[i_1, i_2, ..., i_n] <- value} with one (possibly
          missing) subscript per dimension) only accepts a replacement
          value (a.k.a. right value) that is an array-like object (e.g.
          ordinary array, dgCMatrix object, DelayedArray object, etc...)
          or an ordinary vector (atomic or list) of length 1.

    \item \emph{Linear subassignment} (a.k.a. 1D-style subassignment, that
          is, subassignment of the form \code{x[i] <- value}) only works if
          the subscript \code{i} is a logical DelayedArray object of the same
          dimensions as \code{x} and if the replacement value is an ordinary
          vector (atomic or list) of length 1.

    \item \emph{Filling with a vector}, that is, subassignment of the form
          \code{x[] <- v} where \code{v} is an ordinary vector (atomic or
          list), is only supported if the length of the vector is a divisor
          of \code{nrow(x)}.
  }
  These 3 forms of subassignment are implemented as delayed operations so
  are very light.

  Single value replacement (\code{x[[...]] <- value}) is not supported yet.
}

\seealso{
  \itemize{
    \item \code{\link{showtree}} for \link{DelayedArray} accessors
          \code{nseed}, \code{seed}, and \code{path}.

    \item \code{\link{realize}} for realizing a DelayedArray object in memory
          or on disk.

    \item \code{\link{blockApply}} for more information about block processing
          of array-like objects.

    \item \link{DelayedArray-utils} for common operations on
          DelayedArray objects.

    \item \link{DelayedMatrix-utils} for common operations on
          DelayedMatrix objects.

    \item \link{DelayedArray-stats} for statistical functions on
          DelayedArray objects.

    \item \link{DelayedMatrix-stats} for DelayedMatrix row/col
          summarization.

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
## The seed of a DelayedArray object is **always** treated as a
## "read-only" object so will never be modified by the operations
## we perform on A:
stopifnot(identical(a, seed(A)))
type(A)

## Multi-dimensional single bracket subsetting:
m <- a[11:20 , 5, -3]  # an ordinary matrix
M <- A[11:20 , 5, -3]  # a DelayedMatrix object
stopifnot(identical(m, as.array(M)))

## Linear single bracket subsetting:
A[11:20]
A[A <= 1e-5]
stopifnot(identical(a[a <= 1e-5], A[A <= 1e-5]))

## Subassignment:
A[A < 0.2] <- NA
a[a < 0.2] <- NA
stopifnot(identical(a, as.array(A)))

A[2:5, 1:2, ] <- array(1:40, c(4, 2, 5))
a[2:5, 1:2, ] <- array(1:40, c(4, 2, 5))
stopifnot(identical(a, as.array(A)))

## Other operations:
crazy <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- crazy(a)
head(b)

B <- crazy(A)  # very fast! (all operations are delayed)
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

t(A2)  # transposition is delayed so is very fast and memory-efficient
colSums(A2)

## ---------------------------------------------------------------------
## C. AN HDF5Array OBJECT IS A (PARTICULAR KIND OF) DelayedArray OBJECT
## ---------------------------------------------------------------------
library(HDF5Array)
A3 <- as(a, "HDF5Array")   # write 'a' to an HDF5 file
A3
is(A3, "DelayedArray")     # TRUE
seed(A3)                   # an HDF5ArraySeed object

B3 <- crazy(A3)            # very fast! (all operations are delayed)
B3                         # not an HDF5Array object anymore because
                           # now it carries delayed operations
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
## E. MODIFY THE PATH OF A DelayedArray OBJECT
## ---------------------------------------------------------------------
## This can be useful if the file containing the array data is on a
## shared partition but the exact path to the partition depends on the
## machine from which the data is being accessed.
## For example:

\dontrun{
library(HDF5Array)
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
## F. WRAP A SPARSE MATRIX IN A DelayedArray OBJECT
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

## Murat came up with the following solution that is very fast and
## memory-efficient:
row_normalized_p1 <- Diagonal(x=1/sqrt(Matrix::rowSums(p^2))) %*% p

## With a DelayedArray object, the straight-forward approach uses a
## block processing strategy behind the scene so it doesn't consume
## too much memory.

## First, let's see block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
## and check the automatic block size:
getAutoBlockSize()

row_normalized_P <- P / sqrt(DelayedArray::rowSums(P^2))

## Increasing the block size increases the speed but also memory usage:
setAutoBlockSize(2e8)
row_normalized_P2 <- P / sqrt(DelayedArray::rowSums(P^2))
stopifnot(all.equal(row_normalized_P, row_normalized_P2))

## Back to sparse representation:
DelayedArray:::set_verbose_block_processing(FALSE)
row_normalized_p2 <- as(row_normalized_P, "sparseMatrix")
stopifnot(all.equal(row_normalized_p1, row_normalized_p2))

setAutoBlockSize()
}
}
\keyword{classes}
\keyword{methods}
