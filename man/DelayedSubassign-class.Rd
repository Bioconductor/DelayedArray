\name{DelayedSubassign-class}

\alias{class:DelayedSubassign}
\alias{DelayedSubassign-class}
\alias{DelayedSubassign}

\alias{is_noop,DelayedSubassign-method}
\alias{summary.DelayedSubassign}
\alias{summary,DelayedSubassign-method}

\alias{extract_array,DelayedSubassign-method}

\alias{is_sparse,DelayedSubassign-method}
\alias{extract_sparse_array,DelayedSubassign-method}

\title{DelayedSubassign objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedSubassign class provides a formal representation of a
  \emph{delayed multi-dimensional single bracket subassignment}. It is
  a concrete subclass of the \link{DelayedUnaryIsoOp} virtual class, which
  itself is a subclass of the \link{DelayedUnaryOp} virtual class, which
  itself is a subclass of the \link{DelayedOp} virtual class:
  \preformatted{
                          DelayedOp
                              ^
                              |
                       DelayedUnaryOp
                              ^
                              |
                      DelayedUnaryIsoOp
                              ^
                              |
                      DelayedSubassign
  }

  DelayedSubassign objects are used inside a \link{DelayedArray}
  object to represent the \emph{delayed multi-dimensional single bracket
  subassignments} carried by the object. They're never exposed to the end
  user and are not intended to be manipulated directly.
}

\usage{
\S4method{is_noop}{DelayedSubassign}(x)

\S4method{summary}{DelayedSubassign}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

## DelayedSubassign objects inherit the default dim()
## and dimnames() methods defined for DelayedUnaryIsoOp
## derivatives, but overwite their extract_array() method.

\S4method{extract_array}{DelayedSubassign}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedSubassign}(x)

\S4method{extract_sparse_array}{DelayedSubassign}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedSubassign object.
  }
  \item{index}{
    See \code{?\link[S4Arrays]{extract_array}} in the \pkg{S4Arrays}
    package for a description of the \code{index} argument.
  }
  \item{...}{
    Not used.
  }
}

\seealso{
  \itemize{
    \item \link{DelayedOp} objects.

    \item \code{\link{showtree}} to visualize the nodes and access the
          leaves in the tree of delayed operations carried by a
          \link{DelayedArray} object.

    \item \link[S4Arrays]{extract_array} in the \pkg{S4Arrays} package.

    \item \code{\link{extract_sparse_array}}.
  }
}

\examples{
## DelayedSubassign extends DelayedUnaryIsoOp, which extends
## DelayedUnaryOp, which extends DelayedOp:
extends("DelayedSubassign")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m0 <- matrix(1:30, ncol=5)
M2 <- M1 <- M0 <- DelayedArray(m0)
showtree(M0)

M1[2:5, 5:4] <- 100
showtree(M1)
class(M1@seed)        # a DelayedSubassign object

M2[2:5, 5:4] <- matrix(101:108, ncol=2)
showtree(M2)
class(M2@seed)        # a DelayedSubassign object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------

## DelayedSubassign objects don't propagate sparsity at the moment, that
## is, is_sparse() always returns FALSE on them.

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M1@seed) == "DelayedSubassign")
stopifnot(class(M2@seed) == "DelayedSubassign")
}

\keyword{methods}
