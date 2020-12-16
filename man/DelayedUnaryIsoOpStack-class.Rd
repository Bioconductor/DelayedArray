\name{DelayedUnaryIsoOpStack-class}

\alias{class:DelayedUnaryIsoOpStack}
\alias{DelayedUnaryIsoOpStack-class}
\alias{DelayedUnaryIsoOpStack}

\alias{summary.DelayedUnaryIsoOpStack}
\alias{summary,DelayedUnaryIsoOpStack-method}

\alias{extract_array,DelayedUnaryIsoOpStack-method}

\alias{is_sparse,DelayedUnaryIsoOpStack-method}
\alias{extract_sparse_array,DelayedUnaryIsoOpStack-method}

\title{DelayedUnaryIsoOpStack objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedUnaryIsoOpStack class provides a formal representation of a
  \emph{stack of delayed unary isometric operations}, that is, of a group
  of delayed unary isometric operations stacked (a.k.a. piped) together.
  It is a concrete subclass of the \link{DelayedUnaryIsoOp} virtual class,
  which itself is a subclass of the \link{DelayedUnaryOp} virtual class,
  which itself is a subclass of the \link{DelayedOp} virtual class:
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
                    DelayedUnaryIsoOpStack
  }

  DelayedUnaryIsoOpStack objects are used inside a \link{DelayedArray} object
  to represent groups of delayed unary isometric operations carried by the
  object. They're never exposed to the end user and are not intended to be
  manipulated directly.
}

\usage{
\S4method{summary}{DelayedUnaryIsoOpStack}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

## DelayedUnaryIsoOpStack objects inherit the default dim()
## and dimnames() methods defined for DelayedUnaryIsoOp
## derivatives, but overwite their extract_array() method.

\S4method{extract_array}{DelayedUnaryIsoOpStack}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedUnaryIsoOpStack}(x)

\S4method{extract_sparse_array}{DelayedUnaryIsoOpStack}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedUnaryIsoOpStack object.
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
## DelayedUnaryIsoOpStack extends DelayedUnaryIsoOp, which extends
## DelayedUnaryOp, which extends DelayedOp:
extends("DelayedUnaryIsoOpStack")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m0 <- matrix(runif(12), ncol=3)
M0 <- DelayedArray(m0)
showtree(M0)

M <- log(1 + M0) / 10
showtree(M)
class(M@seed)        # a DelayedUnaryIsoOpStack object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
sm0 <- sparseMatrix(i=c(1, 4), j=c(1, 3), x=c(11, 43), dims=4:3)
SM0 <- DelayedArray(sm0)
showtree(SM0)
is_sparse(SM0)       # TRUE

M1 <- SM0 - 11
showtree(M1)
class(M1@seed)       # a DelayedUnaryIsoOpStack object
is_sparse(M1@seed)   # FALSE

SM2 <- 10 * SM0
showtree(SM2)
class(SM2@seed)      # a DelayedUnaryIsoOpStack object
is_sparse(SM2@seed)  # TRUE

M3 <- SM0 / 0
showtree(M3)
class(M3@seed)       # a DelayedUnaryIsoOpStack object
is_sparse(M3@seed)   # FALSE

SM4 <- log(1 + SM0) / 10
showtree(SM4)
class(SM4@seed)      # a DelayedUnaryIsoOpStack object
is_sparse(SM4@seed)  # TRUE

SM5 <- 2 ^ SM0 - 1
showtree(SM5)
class(SM5@seed)      # a DelayedUnaryIsoOpStack object
is_sparse(SM5@seed)  # TRUE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M@seed) == "DelayedUnaryIsoOpStack")
stopifnot(class(M1@seed) == "DelayedUnaryIsoOpStack")
stopifnot(!is_sparse(M1@seed))
stopifnot(class(SM2@seed) == "DelayedUnaryIsoOpStack")
stopifnot(is_sparse(SM2@seed))
stopifnot(class(M3@seed) == "DelayedUnaryIsoOpStack")
stopifnot(!is_sparse(M3@seed))
stopifnot(class(SM4@seed) == "DelayedUnaryIsoOpStack")
stopifnot(is_sparse(SM4@seed))
stopifnot(class(SM5@seed) == "DelayedUnaryIsoOpStack")
stopifnot(is_sparse(SM5@seed))
}

\keyword{methods}
