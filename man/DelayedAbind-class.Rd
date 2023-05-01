\name{DelayedAbind-class}

\alias{class:DelayedAbind}
\alias{DelayedAbind-class}
\alias{DelayedAbind}

\alias{is_noop,DelayedAbind-method}
\alias{summary.DelayedAbind}
\alias{summary,DelayedAbind-method}

\alias{dim,DelayedAbind-method}
\alias{dimnames,DelayedAbind-method}
\alias{extract_array,DelayedAbind-method}

\alias{is_sparse,DelayedAbind-method}
\alias{extract_sparse_array,DelayedAbind-method}

\alias{updateObject,SeedBinder-method}

\title{DelayedAbind objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedAbind class provides a formal representation of a
  \emph{delayed \code{abind()} operation}. It is a concrete subclass of
  the \link{DelayedNaryOp} virtual class, which itself is a subclass of
  the \link{DelayedOp} virtual class:
  \preformatted{
                          DelayedOp
                              ^
                              |
                        DelayedNaryOp
                              ^
                              |
                        DelayedAbind
  }

  DelayedAbind objects are used inside a \link{DelayedArray} object to
  represent the \emph{delayed \code{abind()} operations} carried by the object.
  They're never exposed to the end user and are not intended to be manipulated
  directly.
}

\usage{
\S4method{is_noop}{DelayedAbind}(x)

\S4method{summary}{DelayedAbind}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{dim}{DelayedAbind}(x)

\S4method{dimnames}{DelayedAbind}(x)

\S4method{extract_array}{DelayedAbind}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedAbind}(x)

\S4method{extract_sparse_array}{DelayedAbind}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedAbind object.
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
## DelayedAbind extends DelayedNaryOp which extends DelayedOp:
extends("DelayedAbind")

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
m1 <- matrix(101:128, ncol=4)
m2 <- matrix(runif(16), ncol=4)
M1 <- DelayedArray(m1)
M2 <- DelayedArray(m2)
showtree(M1)
showtree(M2)

M3 <- rbind(M1, M2)
showtree(M3)
class(M3@seed)        # a DelayedAbind object

M4 <- cbind(t(M1), M2)
showtree(M4)
class(M4@seed)        # a DelayedAbind object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
## DelayedAbind objects always propagate sparsity (granted that all the
## input arrays are sparse).

sm1 <- sparseMatrix(i=c(1, 1, 7, 7), j=c(1, 4, 1, 4),
                    x=c(11, 14, 71, 74), dims=c(7, 4))
SM1 <- DelayedArray(sm1)
sm2 <- sparseMatrix(i=c(1, 1, 4, 4), j=c(1, 4, 1, 4),
                    x=c(11, 14, 41, 44), dims=c(4, 4))
SM2 <- DelayedArray(sm2)
showtree(SM1)
showtree(SM2)
is_sparse(SM1)        # TRUE
is_sparse(SM2)        # TRUE

SM3 <- rbind(SM1, SM2)
showtree(SM3)
class(SM3@seed)       # a DelayedAbind object
is_sparse(SM3@seed)   # TRUE

SM4 <- cbind(SM2, t(SM1))
showtree(SM4)
class(SM4@seed)       # a DelayedAbind object
is_sparse(SM4@seed)   # TRUE

M5 <- rbind(SM2, M1)  # 2nd input array is not sparse!
showtree(M5)
class(M5@seed)        # a DelayedAbind object
is_sparse(M5@seed)    # FALSE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M3@seed) == "DelayedAbind")
stopifnot(class(M4@seed) == "DelayedAbind")
stopifnot(class(SM3@seed) == "DelayedAbind")
stopifnot(is_sparse(SM3@seed))
stopifnot(class(SM4@seed) == "DelayedAbind")
stopifnot(is_sparse(SM4@seed))
stopifnot(class(M5@seed) == "DelayedAbind")
stopifnot(!is_sparse(M5@seed))
}

\keyword{methods}
