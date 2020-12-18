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

\details{
  A DelayedUnaryIsoOpStack object is used to represent the delayed version
  of an operation of the form:
  \preformatted{
    out <- a |> OP1 |> OP2 |> ... |> OPk
  }
  where:
  \itemize{
    \item \code{OP1}, \code{OP2}, ..., \code{OPk} are isometric array
          transformations i.e. operations that return an array with the
          same dimensions as the input array.
    \item \code{a} is the input array.
    \item The output (\code{out}) is an array of same dimensions as \code{a}.
  }

  In addition, each operation (\code{OP}) in the pipe must satisfy the
  property that each value in the output array must be determined **solely**
  by the corresponding value in the input array. In other words:
  \preformatted{
    a |> OP |> `[`(i_1, i_2, ..., i_n)   # i.e. OP(a)[i_1, i_2, ..., i_n]
  }
  must be equal to:
  \preformatted{
    a |> `[`(i_1, i_2, ..., i_n) |> OP   # i.e. OP(a[i_1, i_2, ..., i_n])
  }
  for any valid multidimensional index (i_1, i_2, ..., i_n).

  We refer to this property as the \emph{locality principle}.

  Concrete examples:
  \enumerate{
    \item Things like \code{is.na()}, \code{is.finite()}, logical negation
          (\code{!}), \code{nchar()}, \code{tolower()}.
    \item Most functions in the \link[methods]{Math} and \link[methods]{Math2}
          groups e.g. \code{log()}, \code{sqrt()}, \code{abs()},
          \code{ceiling()}, \code{round()}, etc...
          Notable exceptions are the \code{cum*()} functions (\code{cummin()},
          \code{cummax()}, \code{cumsum()}, and \code{cumprod()}): they don't
          satisfy the \emph{locality principle}.
    \item Operations in the \link[methods]{Ops} group when one operand is
          an array and the other a scalar e.g. \code{a + 10}, \code{2 ^ a},
          \code{a <= 0.5}, etc...
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
