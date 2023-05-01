\name{DelayedAperm-class}

\alias{class:DelayedAperm}
\alias{DelayedAperm-class}
\alias{DelayedAperm}

\alias{is_noop,DelayedAperm-method}
\alias{summary.DelayedAperm}
\alias{summary,DelayedAperm-method}

\alias{dim,DelayedAperm-method}
\alias{dimnames,DelayedAperm-method}
\alias{extract_array,DelayedAperm-method}

\alias{is_sparse,DelayedAperm-method}
\alias{extract_sparse_array,DelayedAperm-method}

\alias{updateObject,SeedDimPicker-method}

\title{DelayedAperm objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedAperm class provides a formal representation of a
  \emph{delayed "extended \code{aperm()}" operation}, that is, of a delayed
  \code{\link[base]{aperm}()} that can drop and/or add \emph{ineffective}
  dimensions. Note that since only \emph{ineffective} dimensions (i.e.
  dimensions with an extent of 1) can be dropped or added, the length of
  the output array is guaranteed to be the same as the length of the input
  array.

  DelayedAperm is a concrete subclass of the \link{DelayedUnaryOp} virtual
  class, which itself is a subclass of the \link{DelayedOp} virtual class:
  \preformatted{
                          DelayedOp
                              ^
                              |
                        DelayedUnaryOp
                              ^
                              |
                         DelayedAperm
  }

  DelayedAperm objects are used inside a \link{DelayedArray} object to
  represent the \emph{delayed "extended \code{aperm()}" operations} carried
  by the object. They're never exposed to the end user and are not intended
  to be manipulated directly.
}

\usage{
\S4method{is_noop}{DelayedAperm}(x)

\S4method{summary}{DelayedAperm}(object, ...)

## ~ ~ ~ Seed contract ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{dim}{DelayedAperm}(x)

\S4method{dimnames}{DelayedAperm}(x)

\S4method{extract_array}{DelayedAperm}(x, index)

## ~ ~ ~ Propagation of sparsity ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

\S4method{is_sparse}{DelayedAperm}(x)

\S4method{extract_sparse_array}{DelayedAperm}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedAperm object.
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
## DelayedAperm extends DelayedUnaryOp which extends DelayedOp:
extends("DelayedAperm")

## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------
a0 <- array(1:20, dim=c(1, 10, 2))
A0 <- DelayedArray(a0)
showtree(A0)

A <- aperm(A0, perm=c(2, 3, 1))
showtree(A)
class(A@seed)       # a DelayedAperm object

M1 <- drop(A0)
showtree(M1)
class(M1@seed)      # a DelayedAperm object

M2 <- t(M1)
showtree(M2)
class(M2@seed)      # a DelayedAperm object

## ---------------------------------------------------------------------
## PROPAGATION OF SPARSITY
## ---------------------------------------------------------------------
## DelayedAperm objects always propagate sparsity.

sa0 <- as(a0, "SparseArraySeed")
SA0 <- DelayedArray(sa0)
showtree(SA0)
is_sparse(SA0)      # TRUE

SA <- aperm(SA0, perm=c(2, 3, 1))
showtree(SA)
class(SA@seed)      # a DelayedAperm object
is_sparse(SA@seed)  # TRUE

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(A@seed) == "DelayedAperm")
stopifnot(class(M1@seed) == "DelayedAperm")
stopifnot(class(M2@seed) == "DelayedAperm")
stopifnot(class(SA@seed) == "DelayedAperm")
stopifnot(is_sparse(SA@seed))
}

\keyword{methods}
