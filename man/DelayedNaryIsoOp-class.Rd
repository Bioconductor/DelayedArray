\name{DelayedNaryIsoOp-class}

\alias{class:DelayedNaryIsoOp}
\alias{DelayedNaryIsoOp-class}
\alias{DelayedNaryIsoOp}

\alias{summary.DelayedNaryIsoOp}
\alias{summary,DelayedNaryIsoOp-method}

\alias{dim,DelayedNaryIsoOp-method}
\alias{dimnames,DelayedNaryIsoOp-method}
\alias{extract_array,DelayedNaryIsoOp-method}

\alias{is_sparse,DelayedNaryIsoOp-method}
\alias{extract_sparse_array,DelayedNaryIsoOp-method}

\alias{updateObject,ConformableSeedCombiner-method}

\title{DelayedNaryIsoOp objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedNaryIsoOp class provides a formal representation of a
  \emph{delayed N-ary isometric operation}. It is a concrete subclass of
  the \link{DelayedNaryOp} virtual class, which itself is a subclass of
  the \link{DelayedOp} virtual class:
  \preformatted{
                          DelayedOp
                              ^
                              |
                        DelayedNaryOp
                              ^
                              |
                       DelayedNaryIsoOp
  }

  DelayedNaryIsoOp objects are used inside a \link{DelayedArray} object to
  represent the \emph{delayed N-ary isometric operation} carried by the object.
  They're never exposed to the end user and are not intended to be manipulated
  directly.
}

\usage{
\S4method{summary}{DelayedNaryIsoOp}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{dim}{DelayedNaryIsoOp}(x)

\S4method{dimnames}{DelayedNaryIsoOp}(x)

\S4method{extract_array}{DelayedNaryIsoOp}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedNaryIsoOp}(x)

\S4method{extract_sparse_array}{DelayedNaryIsoOp}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedNaryIsoOp object.
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
## DelayedNaryIsoOp extends DelayedNaryOp which extends DelayedOp:
extends("DelayedNaryIsoOp")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m1 <- matrix(101:130, ncol=5)
m2 <- matrix(runif(30), ncol=5)
M1 <- DelayedArray(m1)
M2 <- DelayedArray(m2)
showtree(M1)
showtree(M2)

M <- M1 / M2
showtree(M)
class(M@seed)        # a DelayedNaryIsoOp object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
sm1 <- sparseMatrix(i=c(1, 6), j=c(1, 4), x=c(11, 64), dims=6:5)
SM1 <- DelayedArray(sm1)
sm2 <- sparseMatrix(i=c(2, 6), j=c(1, 5), x=c(21, 65), dims=6:5)
SM2 <- DelayedArray(sm2)
showtree(SM1)
showtree(SM2)
is_sparse(SM1)       # TRUE
is_sparse(SM2)       # TRUE

SM3 <- SM1 - SM2
showtree(SM3)
class(SM3@seed)      # a DelayedNaryIsoOp object
is_sparse(SM3@seed)  # TRUE

M4 <- SM1 / SM2
showtree(M4)
class(M4@seed)       # a DelayedNaryIsoOp object
is_sparse(M4@seed)   # FALSE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M@seed) == "DelayedNaryIsoOp")
stopifnot(class(SM3@seed) == "DelayedNaryIsoOp")
stopifnot(is_sparse(SM3@seed))
stopifnot(class(M4@seed) == "DelayedNaryIsoOp")
stopifnot(!is_sparse(M4@seed))
}

\keyword{methods}
