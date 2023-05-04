\name{DelayedSubset-class}

\alias{class:DelayedSubset}
\alias{DelayedSubset-class}
\alias{DelayedSubset}

\alias{is_noop,DelayedSubset-method}
\alias{summary.DelayedSubset}
\alias{summary,DelayedSubset-method}

\alias{dim,DelayedSubset-method}
\alias{dimnames,DelayedSubset-method}
\alias{extract_array,DelayedSubset-method}

\alias{is_sparse,DelayedSubset-method}
\alias{OLD_extract_sparse_array,DelayedSubset-method}

\title{DelayedSubset objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedSubset class provides a formal representation of a \emph{delayed
  multi-dimensional single bracket subsetting operation}. It is a concrete
  subclass of the \link{DelayedUnaryOp} virtual class, which itself is a
  subclass of the \link{DelayedOp} virtual class:
  \preformatted{
                          DelayedOp
                              ^
                              |
                        DelayedUnaryOp
                              ^
                              |
                        DelayedSubset
  }

  DelayedSubset objects are used inside a \link{DelayedArray} object to
  represent the \emph{delayed multi-dimensional single bracket subsetting
  operations} carried by the object. They're never exposed to the end user
  and are not intended to be manipulated directly.
}

\usage{
\S4method{is_noop}{DelayedSubset}(x)

\S4method{summary}{DelayedSubset}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{dim}{DelayedSubset}(x)

\S4method{dimnames}{DelayedSubset}(x)

\S4method{extract_array}{DelayedSubset}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedSubset}(x)

\S4method{OLD_extract_sparse_array}{DelayedSubset}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedSubset object.
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

    \item \code{\link{OLD_extract_sparse_array}}.
  }
}

\examples{
## DelayedSubset extends DelayedUnaryOp which extends DelayedOp:
extends("DelayedSubset")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
a0 <- array(1:60, dim=5:3)
A0 <- DelayedArray(a0)
showtree(A0)

A <- A0[2:1, -4, 3, drop=FALSE]
showtree(A)
class(A@seed)        # a DelayedSubset object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
sm0 <- sparseMatrix(i=c(1, 4), j=c(1, 3), x=c(11, 43), dims=4:3)
SM0 <- DelayedArray(sm0)
showtree(SM0)
is_sparse(SM0)       # TRUE

SM1 <- SM0[-1, 3:2, drop=FALSE]
showtree(SM1)
class(SM1@seed)      # a DelayedSubset object
is_sparse(SM1@seed)  # TRUE

## Duplicated indices break structural sparsity.
M2 <- SM0[-1, c(3:2, 2), drop=FALSE]
showtree(M2)
class(M2@seed)       # a DelayedSubset object
is_sparse(M2@seed)   # FALSE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(A@seed) == "DelayedSubset")
stopifnot(class(SM1@seed) == "DelayedSubset")
stopifnot(is_sparse(SM1@seed))
stopifnot(class(M2@seed) == "DelayedSubset")
stopifnot(!is_sparse(M2@seed))
}

\keyword{methods}
