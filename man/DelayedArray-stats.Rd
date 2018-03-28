\name{DelayedArray-stats}

\alias{DelayedArray-stats}

\alias{dnorm,DelayedArray-method}
\alias{pnorm,DelayedArray-method}
\alias{qnorm,DelayedArray-method}

\alias{dbinom,DelayedArray-method}
\alias{pbinom,DelayedArray-method}
\alias{qbinom,DelayedArray-method}

\alias{dpois,DelayedArray-method}
\alias{ppois,DelayedArray-method}
\alias{qpois,DelayedArray-method}

\alias{dlogis,DelayedArray-method}
\alias{plogis,DelayedArray-method}
\alias{qlogis,DelayedArray-method}

\title{Statistical functions on DelayedArray objects}

\description{
  Statistical functions on \link{DelayedArray} objects.

  All these functions are implemented as delayed operations.
}

\usage{
## --- The Normal Distribution ----- ##

\S4method{dnorm}{DelayedArray}(x, mean=0, sd=1, log=FALSE)
\S4method{pnorm}{DelayedArray}(q, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)
\S4method{qnorm}{DelayedArray}(p, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)

## --- The Binomial Distribution --- ##

\S4method{dbinom}{DelayedArray}(x, size, prob, log=FALSE)
\S4method{pbinom}{DelayedArray}(q, size, prob, lower.tail=TRUE, log.p=FALSE)
\S4method{qbinom}{DelayedArray}(p, size, prob, lower.tail=TRUE, log.p=FALSE)

## --- The Poisson Distribution ---- ##

\S4method{dpois}{DelayedArray}(x, lambda, log=FALSE)
\S4method{ppois}{DelayedArray}(q, lambda, lower.tail=TRUE, log.p=FALSE)
\S4method{qpois}{DelayedArray}(p, lambda, lower.tail=TRUE, log.p=FALSE)

## --- The Logistic Distribution --- ##

\S4method{dlogis}{DelayedArray}(x, location=0, scale=1, log=FALSE)
\S4method{plogis}{DelayedArray}(q, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
\S4method{qlogis}{DelayedArray}(p, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
}

\arguments{
  \item{x, q, p}{
    A \link{DelayedArray} object.
  }
  \item{mean, sd, log, lower.tail, log.p, size, prob, lambda, location, scale}{
    See \code{?stats::\link[stats]{dnorm}}, \code{?stats::\link[stats]{dbinom}},
    \code{?stats::\link[stats]{dpois}}, and \code{?stats::\link[stats]{dlogis}},
    for a description of these arguments.
  }
}

\seealso{
  \itemize{
    \item \code{\link[stats]{dnorm}}, \code{\link[stats]{dbinom}},
          \code{\link[stats]{dpois}}, and \code{\link[stats]{dlogis}}
          in the \pkg{stats} package for the corresponding operations
          on ordinary arrays or matrices.

    \item \link{DelayedArray} objects.

    \item \link[HDF5Array]{HDF5Array} objects in the \pkg{HDF5Array} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(4 * runif(1500000), dim=c(10000, 30, 5))
A <- DelayedArray(a)
A

A2 <- dnorm(A + 1)[ , , -3]  # very fast! (operations are delayed)
A2

a2 <- as.array(A2)           # "realize" 'A2' in memory (as an ordinary
                             # array)

DelayedArray(a2) == A2       # DelayedArray object of type "logical"
stopifnot(all(DelayedArray(a2) == A2))


library(HDF5Array)
A3 <- as(A2, "HDF5Array")    # "realize" 'A2' on disk (as an HDF5Array
                             # object)

A3 == A2                     # DelayedArray object of type "logical"
stopifnot(all(A3 == A2))

## See '?DelayedArray' for general information about DelayedArray objects
## and their "realization".
}
\keyword{methods}
