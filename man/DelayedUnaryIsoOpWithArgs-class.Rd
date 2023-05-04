\name{DelayedUnaryIsoOpWithArgs-class}

\alias{class:DelayedUnaryIsoOpWithArgs}
\alias{DelayedUnaryIsoOpWithArgs-class}
\alias{DelayedUnaryIsoOpWithArgs}

\alias{summary.DelayedUnaryIsoOpWithArgs}
\alias{summary,DelayedUnaryIsoOpWithArgs-method}

\alias{extract_array,DelayedUnaryIsoOpWithArgs-method}

\alias{is_sparse,DelayedUnaryIsoOpWithArgs-method}
\alias{OLD_extract_sparse_array,DelayedUnaryIsoOpWithArgs-method}

\title{DelayedUnaryIsoOpWithArgs objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  The DelayedUnaryIsoOpWithArgs class provides a formal representation of a
  \emph{delayed unary isometric operation with vector-like arguments going
  along the dimensions of the input array}. It is a concrete subclass of
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
                  DelayedUnaryIsoOpWithArgs
  }

  DelayedUnaryIsoOpWithArgs objects are used inside a \link{DelayedArray}
  object to represent the \emph{delayed unary isometric operations with
  vector-like arguments going along the dimensions of the input array} carried
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

\S4method{OLD_extract_sparse_array}{DelayedUnaryIsoOpWithArgs}(x, index)
}

\arguments{
  \item{x, object}{
    A DelayedUnaryIsoOpWithArgs object.
  }
  \item{index}{
    See \code{?\link[S4Arrays]{extract_array}} in the \pkg{S4Arrays}
    package for a description of the \code{index} argument.
  }
  \item{...}{
    Not used.
  }
}

\details{
  A DelayedUnaryIsoOpWithArgs object is used to represent the delayed version
  of an operation of the form:
  \preformatted{
    out <- OP(L1, L2, ..., a, R1, R2, ...)
  }
  where:
  \itemize{
    \item \code{OP} is an isometric array transformation i.e. a transformation
          that returns an array with the same dimensions as the input array.
    \item \code{a} is the input array.
    \item \code{L1}, \code{L2}, etc... are the left arguments.
    \item \code{R1}, \code{R2}, etc... are the right arguments.
    \item The output (\code{out}) is an array of same dimensions as \code{a}.
  }
  Some of the arguments (left or right) can go along the dimensions of the
  input array. For example if \code{a} is a 12 x 150 x 5 array, argument
  \code{L2} is considered to go along the 3rd dimension if its length is 5
  and if the result of:
  \preformatted{
    OP(L1, L2[k], ..., a[ , , k, drop=FALSE], R1, R2, ...)
  }
  is the same as \code{out[ , , k, drop=FALSE]} for any index \code{k}.

  More generally speaking, if, say, arguments \code{L2}, \code{L3}, \code{R1},
  and \code{R2} go along the 3rd, 1st, 2nd, and 1st dimensions, respectively,
  then each value in the output array (\code{a[i, j, k]}) must be determined
  solely by the corresponding values in the input array (\code{a[i, j, k]})
  and arguments (\code{L2[k]}, \code{L3[i]}, \code{R1[j]}, \code{R2[i]}).
  In other words, \code{out[i, j, k]} must be equal to:
  \preformatted{
    OP(L1, L2[k], L3[i], ..., a[i, j, k], R1[j], R2[i], ...)
  }
  for any 1 <= \code{i} <= 12, 1 <= \code{j} <= 150, and 1 <= \code{k} <= 5.

  We refer to this property as the \emph{locality principle}.

  Concrete examples:
  \enumerate{
    \item Addition (or any operation in the \link[methods]{Ops} group)
          of an array \code{a} and an atomic vector \code{v} of length
          \code{dim(a)[[1]]}:
          \itemize{
            \item \code{`+`(a, v)}: \code{OP} is \code{`+`}, right argument
                  goes along the 1st dimension.
            \item \code{`<=`(a, v)}: \code{OP} is \code{`<=`}, right argument
                  goes along the 1st dimension.
            \item \code{`&`(v, a)}: \code{OP} is \code{`&`}, left argument
                  goes along the 1st dimension.
          }
    \item \code{scale(x, center=v1, scale=v2)}: \code{OP} is \code{scale},
          right arguments \code{center} and \code{scale} go along the 2nd
          dimension.
  }

  Note that if \code{OP} has no argument that goes along a dimension of
  the input array, then the delayed operation is better represented with
  a \link{DelayedUnaryIsoOpStack} object.
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

SM2 <- SM0 * 101:104
showtree(SM2)
class(SM2@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(SM2@seed)  # TRUE

SM3 <- SM0 * c(101:103, 0)
showtree(SM3)
class(SM3@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(SM3@seed)  # TRUE

M4 <- SM0 * c(101:103, NA)
showtree(M4)
class(M4@seed)       # a DelayedUnaryIsoOpWithArgs object
is_sparse(M4@seed)   # FALSE

M5 <- SM0 * c(101:103, Inf)
showtree(M5)
class(M5@seed)       # a DelayedUnaryIsoOpWithArgs object
is_sparse(M5@seed)   # FALSE

SM6 <- SM0 / 101:104
showtree(SM6)
class(SM6@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(SM6@seed)  # TRUE

M7 <- SM0 / c(101:103, 0)
showtree(M7)
class(M7@seed)       # a DelayedUnaryIsoOpWithArgs object
is_sparse(M7@seed)   # FALSE

M8 <- SM0 / c(101:103, NA)
showtree(M8)
class(M8@seed)       # a DelayedUnaryIsoOpWithArgs object
is_sparse(M8@seed)   # FALSE

SM9 <- SM0 / c(101:103, Inf)
showtree(SM9)
class(SM9@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(SM9@seed)  # TRUE

M10 <- 101:104 / SM0
showtree(M10)
class(M10@seed)      # a DelayedUnaryIsoOpWithArgs object
is_sparse(M10@seed)  # FALSE

## ---------------------------------------------------------------------
## ADVANCED EXAMPLE
## ---------------------------------------------------------------------
## Not ready yet!
#op <- DelayedArray:::new_DelayedUnaryIsoOpWithArgs(m0,
#          scale,
#          Rargs=list(center=c(1, 0, 100), scale=c(10, 1, 1)),
#          Ralong=c(2, 2))

## ---------------------------------------------------------------------
## SANITY CHECKS
## ---------------------------------------------------------------------
stopifnot(class(M@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(class(M1@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M1@seed))
stopifnot(class(SM2@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(is_sparse(SM2@seed))
stopifnot(class(SM3@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(is_sparse(SM3@seed))
stopifnot(class(M4@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M4@seed))
stopifnot(class(M5@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M5@seed))
stopifnot(class(SM6@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(is_sparse(SM6@seed))
stopifnot(class(M7@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M7@seed))
stopifnot(class(M8@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M8@seed))
stopifnot(class(SM9@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(is_sparse(SM9@seed))
stopifnot(class(M10@seed) == "DelayedUnaryIsoOpWithArgs")
stopifnot(!is_sparse(M10@seed))
}

\keyword{methods}
