\name{DelayedUnaryIsoOpWithArgs-class}

\alias{class:DelayedUnaryIsoOpWithArgs}
\alias{DelayedUnaryIsoOpWithArgs-class}
\alias{DelayedUnaryIsoOpWithArgs}

\alias{summary.DelayedUnaryIsoOpWithArgs}
\alias{summary,DelayedUnaryIsoOpWithArgs-method}

\alias{extract_array,DelayedUnaryIsoOpWithArgs-method}

\alias{is_sparse,DelayedUnaryIsoOpWithArgs-method}
\alias{extract_sparse_array,DelayedUnaryIsoOpWithArgs-method}

\title{DelayedUnaryIsoOpWithArgs objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedUnaryIsoOpWithArgs class provides a formal representation of a
  \emph{delayed unary isometric operation with vector-like arguments along
  the dimensions of the input array}. It is a concrete subclass of the
  \link{DelayedUnaryIsoOp} virtual class, which itself is a subclass of
  the \link{DelayedUnaryOp} virtual class, which itself is a subclass of
  the \link{DelayedOp} virtual class:
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
                  DelayedUnaryIsoOpWithArgs
  }

  DelayedUnaryIsoOpWithArgs objects are used inside a \link{DelayedArray}
  object to represent the \emph{delayed unary isometric operations with
  vector-like arguments along the dimensions of the input array} carried
  by the object. They're never exposed to the end user and are not intended
  to be manipulated directly.
}

\usage{
\S4method{summary}{DelayedUnaryIsoOpWithArgs}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

## DelayedUnaryIsoOpWithArgs objects inherit the default dim()
## and dimnames() methods defined for DelayedUnaryIsoOp
## derivatives, but overwite their extract_array() method.

\S4method{extract_array}{DelayedUnaryIsoOpWithArgs}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedUnaryIsoOpWithArgs}(x)

\S4method{extract_sparse_array}{DelayedUnaryIsoOpWithArgs}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedUnaryIsoOpWithArgs object.
  }
  \item{index}{
    See \code{?\link{extract_array}} for a description of the \code{index}
    argument.
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

    \item \code{\link{extract_array}} and \code{\link{extract_sparse_array}}.
  }
}

\examples{
## DelayedUnaryIsoOpWithArgs extends DelayedUnaryIsoOp, which extends
## DelayedUnaryOp, which extends DelayedOp:
extends("DelayedUnaryIsoOpWithArgs")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m0 <- matrix(runif(12), ncol=3)
M0 <- DelayedArray(m0)
showtree(M0)

M <- M0 + 101:104
showtree(M)
class(M@seed)        # a DelayedUnaryIsoOpWithArgs object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
sm0 <- sparseMatrix(i=c(1, 4), j=c(1, 3), x=c(11, 43), dims=4:3)
SM0 <- DelayedArray(sm0)
showtree(SM0)
is_sparse(SM0)       # TRUE

M1 <- SM0 + 101:104
showtree(M1)
class(M1@seed)       # a DelayedUnaryIsoOpWithArgs object
is_sparse(M1@seed)   # FALSE

SM2 <- SM0 / 101:104
showtree(SM2)
class(SM2@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(SM2@seed)  # TRUE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(class(M1@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M1@seed))
stopifnot(class(SM2@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(is_sparse(SM2@seed))
}

\keyword{methods}
