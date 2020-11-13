\name{DelayedArray-utils}

\alias{DelayedArray-utils}

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

\alias{is.na,DelayedArray-method}
\alias{is.finite,DelayedArray-method}
\alias{is.infinite,DelayedArray-method}
\alias{is.nan,DelayedArray-method}
\alias{!,DelayedArray-method}

\alias{type<-,DelayedArray,ANY-method}

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

\alias{unique}
\alias{unique.DelayedArray}
\alias{unique,DelayedArray-method}
\alias{table}
\alias{table,DelayedArray-method}

\alias{range}
\alias{range.DelayedArray}
\alias{range,DelayedArray-method}

\alias{mean}
\alias{mean.DelayedArray}
\alias{mean,DelayedArray-method}

\alias{apply}
\alias{apply,DelayedArray-method}

\alias{sweep}
\alias{sweep,DelayedArray-method}

\alias{scale}
\alias{scale.DelayedMatrix}
\alias{scale,DelayedMatrix-method}

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
    \item \code{!}
    \item \code{is.na}, \code{is.finite}, \code{is.infinite}, \code{is.nan}
    \item \code{type<-}
    \item \code{lengths}
    \item \code{nchar}, \code{tolower}, \code{toupper},
          \code{grepl}, \code{sub}, \code{gsub}
    \item \code{pmax2} and \code{pmin2}
    \item \code{\link[base]{sweep}}
    \item \code{\link[base]{scale}} (when the supplied \code{center} and
          \code{scale} are not TRUE)
    \item statistical functions like \code{dnorm}, \code{dbinom}, \code{dpois},
          and \code{dlogis} (for the Normal, Binomial, Poisson, and Logistic
          distribution, respectively) and related functions (documented in
          \link{DelayedArray-stats})
  }

  Block-processed operations:
  \itemize{
    \item \code{anyNA}, \code{which}
    \item \code{unique}, \code{table}
    \item all the members of the \code{\link[methods]{Summary}} group
    \item \code{mean}
    \item \code{apply}
  }

  Mix delayed and block-processed operations:
  \itemize{
    \item \code{\link[base]{scale}} (when the supplied \code{center} and/or
          \code{scale} are TRUE)
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
          \code{\link[base]{table}}, \code{\link[base]{mean}},
          \code{\link[base]{apply}}, and \code{\link[base]{\%*\%}} in the
          \pkg{base} package for the corresponding operations on ordinary
          arrays or matrices.

    \item \link{DelayedMatrix-utils} for common operations on
          \link{DelayedMatrix} objects.

    \item \link{DelayedArray-stats} for statistical functions on
          \link{DelayedArray} objects.

    \item \link{DelayedMatrix-stats} for \link{DelayedMatrix} row/col
          summarization.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package
          for the members of the \code{\link[methods]{Ops}},
          \code{\link[methods]{Math}}, and \code{\link[methods]{Math2}} groups.

    \item \code{\link[base]{sweep}} and \code{\link[base]{scale}} in the
          \pkg{base} package.
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

M12 <- rbind(M1, t(M2))            # delayed
M12
colMeans(M12)                      # block-processed

## On objects with more than 2 dimensions:
example(arbind)  # to create arrays a1, a2, a3

A1 <- DelayedArray(a1)
A2 <- DelayedArray(a2)
A3 <- DelayedArray(a3)
A123 <- rbind(A1, A2, A3)          # delayed
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

M1 >= 0.5 & M1 < 0.75              # delayed
log(M1)                            # delayed
pmax2(M2, 0)                       # delayed

type(M2) <- "integer"              # delayed
M2

## table() is block-processed:
a4 <- array(sample(50L, 2000000L, replace=TRUE), c(200, 4, 2500))
A4 <- as(a4, "HDF5Array")
table(A4)
a5 <- array(sample(20L, 2000000L, replace=TRUE), c(200, 4, 2500))
A5 <- as(a5, "HDF5Array")
table(A5)

A4 - 2 * A5                        # delayed
table(A4 - 2 * A5)                 # block-processed

## range() is block-processed:
range(A4 - 2 * A5)
range(M1)

cmeans <- colMeans(M2)             # block-processed
sweep(M2, 2, cmeans)               # delayed
scale(M2)                          # delayed & block-processed
scale(M2, center=FALSE, scale=10)  # delayed
}
\keyword{methods}
