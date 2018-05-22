\name{DelayedArray-utils}

\alias{DelayedArray-utils}

% DelayedArray utils

\alias{arbind,DelayedArray-method}
\alias{acbind,DelayedArray-method}
\alias{rbind}
\alias{cbind}
\alias{rbind,DelayedArray-method}
\alias{cbind,DelayedArray-method}
\alias{bindROWS,DelayedArray-method}

\alias{+,DelayedArray,missing-method}
\alias{-,DelayedArray,missing-method}

\alias{pmax2}
\alias{pmin2}
\alias{pmax2,ANY,ANY-method}
\alias{pmin2,ANY,ANY-method}
\alias{pmax2,DelayedArray,vector-method}
\alias{pmin2,DelayedArray,vector-method}
\alias{pmax2,vector,DelayedArray-method}
\alias{pmin2,vector,DelayedArray-method}
\alias{pmax2,DelayedArray,DelayedArray-method}
\alias{pmin2,DelayedArray,DelayedArray-method}

\alias{sweep,DelayedArray-method}

\alias{is.na,DelayedArray-method}
\alias{is.finite,DelayedArray-method}
\alias{is.infinite,DelayedArray-method}
\alias{is.nan,DelayedArray-method}
\alias{!,DelayedArray-method}

\alias{lengths,DelayedArray-method}

\alias{nchar,DelayedArray-method}
\alias{tolower,DelayedArray-method}
\alias{toupper,DelayedArray-method}

\alias{round,DelayedArray-method}
\alias{signif,DelayedArray-method}

\alias{grepl,ANY,DelayedArray-method}
\alias{sub,ANY,ANY,DelayedArray-method}
\alias{gsub,ANY,ANY,DelayedArray-method}

\alias{anyNA,DelayedArray-method}
\alias{which,DelayedArray-method}

\alias{mean.DelayedArray}
\alias{mean,DelayedArray-method}

\alias{apply}
\alias{apply,DelayedArray-method}

% DelayedMatrix utils

\alias{t.DelayedMatrix}
\alias{t,DelayedMatrix-method}

\alias{\%*\%}
\alias{\%*\%,DelayedMatrix,matrix-method}
\alias{\%*\%,matrix,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}

% DelayedMatrix row/col summarization

\alias{rowSums}
\alias{rowSums,DelayedMatrix-method}
\alias{colSums}
\alias{colSums,DelayedMatrix-method}

\alias{rowMeans}
\alias{rowMeans,DelayedMatrix-method}
\alias{colMeans}
\alias{colMeans,DelayedMatrix-method}

\alias{rowMaxs}
\alias{rowMaxs,DelayedMatrix-method}
\alias{colMaxs}
\alias{colMaxs,DelayedMatrix-method}

\alias{rowMins}
\alias{rowMins,DelayedMatrix-method}
\alias{colMins}
\alias{colMins,DelayedMatrix-method}

\alias{rowRanges}
\alias{rowRanges,DelayedMatrix-method}
\alias{colRanges}
\alias{colRanges,DelayedMatrix-method}

\title{Common operations on DelayedArray objects}

\description{
  Common operations on \link{DelayedArray} objects.
}

\details{
  The operations currently supported on \link{DelayedArray} objects are:

  Delayed operations:
  \itemize{
    \item \code{rbind} and \code{cbind}
    \item all the members of the \code{\link[methods]{Ops}},
          \code{\link[methods]{Math}}, and \code{\link[methods]{Math2}} groups
    \item \code{sweep}
    \item \code{!}
    \item \code{is.na}, \code{is.finite}, \code{is.infinite}, \code{is.nan}
    \item \code{lengths}
    \item \code{nchar}, \code{tolower}, \code{toupper},
          \code{grepl}, \code{sub}, \code{gsub}
    \item \code{pmax2} and \code{pmin2}
    \item \code{t}
    \item statistical functions like \code{dnorm}, \code{dbinom}, \code{dpois},
          and \code{dlogis} (for the Normal, Binomial, Poisson, and Logistic
          distribution, respectively) and related functions (documented in
          \link{DelayedArray-stats})
  }

  Block-processed operations:
  \itemize{
    \item \code{anyNA}, \code{which}
    \item all the members of the \code{\link[methods]{Summary}} group
    \item \code{mean}
    \item \code{apply}
    \item matrix multiplication (\%*\%) of an ordinary matrix by a
          \link{DelayedMatrix} object
    \item matrix row/col summarization [\link{DelayedMatrix} objects only]:
          \code{rowSums}, \code{colSums}, \code{rowMeans}, \code{colMeans},
          \code{rowMaxs}, \code{colMaxs}, \code{rowMins}, \code{colMins},
          \code{rowRanges}, and \code{colRanges}
  }
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in the \pkg{base} package for
          rbind/cbind'ing ordinary arrays.

    \item \code{\link{arbind}} and \code{\link{acbind}} in this package
          (\pkg{DelayedArray}) for binding ordinary arrays of arbitrary
          dimensions along their rows or columns.

    \item \code{\link[base]{is.na}}, \code{\link[base]{!}},
          \code{\link[base]{mean}}, \code{\link[base]{apply}},
          and \code{\link[base]{\%*\%}} in the \pkg{base} package for the
          corresponding operations on ordinary arrays or matrices.

    \item \code{\link[base]{rowSums}} in the \pkg{base} package and
          \code{\link[matrixStats]{rowMaxs}} in the \pkg{matrixStats} package
          for row/col summarization of an ordinary matrix.

    \item \code{\link{setRealizationBackend}} for how to set a
          \emph{realization backend}.

    \item \code{\link[HDF5Array]{writeHDF5Array}} in the \pkg{HDF5Array}
          package for writing an array-like object to an HDF5 file and other
          low-level utilities to control the location of automatically created
          HDF5 datasets.

    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-stats} for statistical functions on
          DelayedArray objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package
          for the members of the \code{\link[methods]{Ops}},
          \code{\link[methods]{Math}}, and \code{\link[methods]{Math2}} groups.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## BIND DelayedArray OBJECTS
## ---------------------------------------------------------------------
## DelayedArray objects can be bound along their 1st (rows) or 2nd
## (columns) dimension with rbind() or cbind(). These operations are
## equivalent to arbind() and acbind(), respectively, and are all
## delayed.

## On 2D objects:
library(HDF5Array)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M12 <- rbind(M1, t(M2))        # delayed
M12
colMeans(M12)                  # block-processed

## On objects with more than 2 dimensions:
example(arbind)  # to create arrays a1, a2, a3

A1 <- DelayedArray(a1)
A2 <- DelayedArray(a2)
A3 <- DelayedArray(a3)
A123 <- rbind(A1, A2, A3)      # delayed
A123

## On 1D objects:
v1 <- array(11:15, 5, dimnames=list(LETTERS[1:5]))
v2 <- array(letters[1:3])
V1 <- DelayedArray(v1)
V2 <- DelayedArray(v2)
V12 <- rbind(V1, V2)
V12

\dontrun{cbind(V1, V2)  # Error! (the objects to cbind() must have at least 2
               # dimensions)
}

## Note that base::rbind() and base::cbind() do something completely
## different on ordinary arrays that are not matrices. They treat them
## as if they were vectors:
rbind(a1, a2, a3)
cbind(a1, a2, a3)
rbind(v1, v2)
cbind(v1, v2)

## Also note that DelayedArray objects of arbitrary dimensions can be
## stored inside a DataFrame object as long as they all have the same
## first dimension (nrow()):
DF <- DataFrame(M=I(tail(M1, n=5)), A=I(A3), V=I(V1))
DF[-3, ]
DF2 <- rbind(DF, DF)
DF2$V

## Sanity checks:
m1 <- as.matrix(M1)
m2 <- as.matrix(M2)
stopifnot(identical(rbind(m1, t(m2)), as.matrix(M12)))
stopifnot(identical(arbind(a1, a2, a3), as.array(A123)))
stopifnot(identical(arbind(v1, v2), as.array(V12)))
stopifnot(identical(rbind(DF$M, DF$M), DF2$M))
stopifnot(identical(rbind(DF$A, DF$A), DF2$A))
stopifnot(identical(rbind(DF$V, DF$V), DF2$V))

## ---------------------------------------------------------------------
## MORE OPERATIONS
## ---------------------------------------------------------------------

range(M1)                      # block-processed
M1 >= 0.5 & M1 < 0.75          # delayed
log(M1)                        # delayed

pmax2(M2, 0)                   # delayed

sweep(M2, 2, colMeans(M2))     # delayed

## ---------------------------------------------------------------------
## MATRIX MULTIPLICATION
## ---------------------------------------------------------------------

## Matrix multiplication is not delayed: the output matrix is realized
## block by block. The current "realization backend" controls where
## realization happens e.g. in memory if set to NULL or in an HDF5 file
## if set to "HDF5Array". See '?realize' for more information about
## "realization backends".
## The output matrix is returned as a DelayedMatrix object with no delayed
## operations on it. The exact class of the object depends on the backend
## e.g. it will be HDF5Matrix with "HDF5Array" backend.

m <- matrix(runif(50000), ncol=nrow(M1))

## Set backend to NULL for in-memory realization:
setRealizationBackend()
P1 <- m \%*\% M1
P1

## Set backend to HDF5Array for realization in HDF5 file:
setRealizationBackend("HDF5Array")

## With the HDF5Array backend, the output matrix will be written to an
## automatic location on disk:
getHDF5DumpFile()  # HDF5 file where the output matrix will be written
lsHDF5DumpFile()

P2 <- m \%*\% M1
P2

lsHDF5DumpFile()

## Use setHDF5DumpFile() and setHDF5DumpName() from the HDF5Array package
## to control the location of automatically created HDF5 datasets.

stopifnot(identical(as.array(P1), as.array(P2)))

## ---------------------------------------------------------------------
## MATRIX ROW/COL SUMMARIZATION
## ---------------------------------------------------------------------

## All these operations are block-processed.

rowSums(M1)
colSums(M1)

rowMeans(M1)
colMeans(M1)

rmaxs <- rowMaxs(M1)
cmaxs <- colMaxs(M1)

rmins <- rowMins(M1)
cmins <- colMins(M1)

rranges <- rowRanges(M1)
cranges <- colRanges(M1)

stopifnot(identical(cbind(rmins, rmaxs, deparse.level=0), rranges))
stopifnot(identical(cbind(cmins, cmaxs, deparse.level=0), cranges))
}
\keyword{methods}
