\name{DelayedSetDimnames-class}

\alias{class:DelayedSetDimnames}
\alias{DelayedSetDimnames-class}
\alias{DelayedSetDimnames}

\alias{is_noop,DelayedSetDimnames-method}
\alias{summary.DelayedSetDimnames}
\alias{summary,DelayedSetDimnames-method}

\alias{dimnames,DelayedSetDimnames-method}

\alias{updateObject,DelayedDimnames-method}

\title{DelayedSetDimnames objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedSetDimnames class provides a formal representation of a
  \emph{delayed "set dimnames" operation}. It is a concrete subclass of
  the \link{DelayedUnaryIsoOp} virtual class, which itself is a subclass of
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
                      DelayedSetDimnames
  }

  DelayedSetDimnames objects are used inside a \link{DelayedArray}
  object to represent the \emph{delayed "set dimnames" operations}
  carried by the object. They're never exposed to the end user and
  are not intended to be manipulated directly.
}

\usage{
\S4method{is_noop}{DelayedSetDimnames}(x)

\S4method{summary}{DelayedSetDimnames}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

## DelayedSetDimnames objects inherit the default dim()
## and extract_array() methods defined for DelayedUnaryIsoOp
## derivatives, but overwite their dimnames() method.

\S4method{dimnames}{DelayedSetDimnames}(x)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

## DelayedSetDimnames objects inherit the default
## is_sparse() and extract_sparse_array() methods defined
## for DelayedUnaryIsoOp derivatives.
}

\arguments{
  \item{x, object}{
    A DelayedSetDimnames object.
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
  }
}

\examples{
## DelayedSetDimnames extends DelayedUnaryIsoOp, which extends
## DelayedUnaryOp, which extends DelayedOp:
extends("DelayedSetDimnames")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m0 <- matrix(1:30, ncol=5, dimnames=list(letters[1:6], NULL))
M2 <- M1 <- M0 <- DelayedArray(m0)
showtree(M0)

dimnames(M1) <- list(NULL, LETTERS[1:5])
showtree(M1)
class(M1@seed)      # a DelayedSetDimnames object

colnames(M2) <- LETTERS[1:5]
showtree(M2)
class(M2@seed)      # a DelayedSetDimnames object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
## DelayedSetDimnames objects always propagate sparsity.

sm0 <- sparseMatrix(i=c(1, 4), j=c(1, 3), x=c(11, 43), dims=4:3)
SM <- SM0 <- DelayedArray(sm0)
showtree(SM0)
is_sparse(SM0)      # TRUE

dimnames(SM) <- list(letters[1:4], LETTERS[1:3])
showtree(SM)
class(SM@seed)      # a DelayedSetDimnames object
is_sparse(SM@seed)  # TRUE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M1@seed) == "DelayedSetDimnames")
stopifnot(class(M2@seed) == "DelayedSetDimnames")
stopifnot(class(SM@seed) == "DelayedSetDimnames")
stopifnot(is_sparse(SM@seed))
}

\keyword{methods}
