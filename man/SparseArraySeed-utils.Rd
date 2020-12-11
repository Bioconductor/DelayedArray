\name{SparseArraySeed-utils}

\alias{SparseArraySeed-utils}

\alias{rbind,SparseArraySeed-method}
\alias{cbind,SparseArraySeed-method}

\alias{is.na,SparseArraySeed-method}
\alias{is.infinite,SparseArraySeed-method}
\alias{is.nan,SparseArraySeed-method}
\alias{tolower,SparseArraySeed-method}
\alias{toupper,SparseArraySeed-method}
\alias{nchar,SparseArraySeed-method}

\alias{anyNA,SparseArraySeed-method}
\alias{which,SparseArraySeed-method}

\alias{max,SparseArraySeed-method}
\alias{min,SparseArraySeed-method}
\alias{range.SparseArraySeed}
\alias{range,SparseArraySeed-method}
\alias{sum,SparseArraySeed-method}
\alias{prod,SparseArraySeed-method}
\alias{any,SparseArraySeed-method}
\alias{all,SparseArraySeed-method}

\alias{mean.SparseArraySeed}
\alias{mean,SparseArraySeed-method}

\title{Operate natively on SparseArraySeed objects}

\description{
   Some utilities to operate natively on \link{SparseArraySeed} objects.
   Mostly for internal use by the \pkg{DelayedArray} package e.g. they
   support block processed methods for sparse \link{DelayedArray} objects
   like \code{sum()}, \code{mean()}, \code{which()}, etc...
}

\usage{
\S4method{is.na}{SparseArraySeed}(x)

\S4method{is.infinite}{SparseArraySeed}(x)

\S4method{is.nan}{SparseArraySeed}(x)

\S4method{tolower}{SparseArraySeed}(x)

\S4method{toupper}{SparseArraySeed}(x)

\S4method{nchar}{SparseArraySeed}(x, type="chars", allowNA=FALSE, keepNA=NA)

\S4method{anyNA}{SparseArraySeed}(x, recursive=FALSE)

\S4method{which}{SparseArraySeed}(x, arr.ind=FALSE, useNames=TRUE)


## <>-<>-<> "Summary" group generic <>-<>-<>

\S4method{max}{SparseArraySeed}(x, ..., na.rm=FALSE)

\S4method{min}{SparseArraySeed}(x, ..., na.rm=FALSE)

\S4method{range}{SparseArraySeed}(x, ..., finite=FALSE, na.rm=FALSE)

\S4method{sum}{SparseArraySeed}(x, ..., na.rm=FALSE)

\S4method{prod}{SparseArraySeed}(x, ..., na.rm=FALSE)

\S4method{any}{SparseArraySeed}(x, ..., na.rm=FALSE)

\S4method{all}{SparseArraySeed}(x, ..., na.rm=FALSE)


## <>-<>-<> others <>-<>-<>

\S4method{mean}{SparseArraySeed}(x, na.rm=FALSE)
}

\arguments{
  \item{x}{
    A \link{SparseArraySeed} object.
  }
  \item{type, allowNA, keepNA}{
    See \code{?base::\link[base]{nchar}} for a description of
    these arguments.
  }
  \item{recursive}{
    See \code{?base::\link[base]{anyNA}} for a description of
    this argument.
  }
  \item{arr.ind}{
    See \code{?base::\link[base]{which}} for a description of
    this argument.
  }
  \item{useNames}{
    Ignored.
  }
  \item{...}{
    Unsupported.
  }
  \item{na.rm}{
    \code{TRUE} or \code{FALSE} (the default).
    Should \code{NA}'s and \code{NaN}'s be removed?
  }
  \item{finite}{
    \code{TRUE} or \code{FALSE} (the default).
    Should non-finite values be removed?
  }
}

\value{
  See corresponding functions in the \pkg{base} package.
}

\seealso{
  \itemize{
    \item \link{SparseArraySeed} objects.
  }
}

\examples{
## Create a SparseArraySeed object:
dim1 <- 5:3
nzindex1 <- Lindex2Mindex(sample(60, 14), 5:3)
sas1 <- SparseArraySeed(dim1, nzindex1, nzdata=sample(0:13))

## Apply native operations:
sum(sas1)
range(sas1)
mean(sas1)

## Sanity checks:
stopifnot(identical(sum(as.array(sas1)), sum(sas1)))
stopifnot(identical(range(as.array(sas1)), range(sas1)))
stopifnot(identical(mean(as.array(sas1)), mean(sas1)))
}
\keyword{methods}
