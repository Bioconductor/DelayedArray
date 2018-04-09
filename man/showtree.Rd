\name{showtree}

% showtree()
\alias{showtree}
\alias{show,DelayedOp-method}

% simplify()
\alias{simplify}
\alias{simplify,DelayedSubset-method}
\alias{simplify,DelayedAperm-method}
\alias{simplify,DelayedUnaryIsoOp-method}
\alias{simplify,DelayedDimnames-method}
\alias{simplify,DelayedArray-method}

% summarizeMappingToSeed()
\alias{summarizeMappingToSeed}
\alias{summarizeMappingToSeed,ANY-method}
\alias{summarizeMappingToSeed,DelayedArray-method}

\title{Visualize and simplify a tree of delayed operations}

\description{
  NOTE: The tools documented in this man page are primarily intended for
  developers. End users of \link{DelayedArray} objects will typically not
  need them.

  \code{showtree} can be used to visualize the tree of delayed operations
  carried by a \link{DelayedArray} object.

  \code{simplify} can be used to simplify this tree.

  \code{summarizeMappingToSeed} summarizes the mapping between the elements
  of a \link{DelayedArray} object and the elements of its seed.
}

\usage{
showtree(x, show.node.dim=TRUE)
simplify(x, incremental=FALSE)
summarizeMappingToSeed(x, as.DelayedOp=FALSE)
}

\arguments{
  \item{x}{
    Typically a \link{DelayedArray} object but can also be a \link{DelayedOp}
    object.

    Additionally \code{showtree} accepts a list where each element is a
    \link{DelayedArray} or \link{DelayedOp} object.
  }
  \item{show.node.dim}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} (the default), the nodes
    dimensions and data type are displayed.
  }
  \item{incremental}{
    For internal use.
  }
  \item{as.DelayedOp}{
    \code{FALSE} by default. If \code{TRUE}, \code{summarizeMappingToSeed}
    returns a linear tree with 2 \link{DelayedOp} nodes.
    See details below.
  }
}

\details{
  \code{summarizeMappingToSeed} is only supported on a \link{DelayedArray}
  object \code{x} with a single seed i.e. if \code{nseed(x) == 1}.

  The mapping between the elements of \code{x} and the elements of its seed
  is affected by the following delayed operations carried by \code{x}:
  \code{[}, \code{drop()}, and \code{aperm()}.
  \code{x} can carry any number of each of these operations in any order but
  their net result can be described by a \emph{net subsetting} followed by
  a \emph{net dimension rearrangement}.

  \code{summarizeMappingToSeed(x)} returns an object that represents the
  \emph{net subsetting} and \emph{net dimension rearrangement}. This object
  can be used to map the elements of \code{x} to their corresponding element
  in \code{seed(x)}.

  By default (i.e. if \code{as.DelayedOp} is \code{FALSE}), this object is
  a list of subscripts that describes the \emph{net subsetting}. The list
  contains one subscript per dimension in the seed. Each subscript can
  be either a vector of positive integers or a \code{NULL}. A \code{NULL}
  indicates a missing subscript.
  The \emph{net dimension rearrangement} is described in a \code{dimmap}
  attribute added to the list. This attribute is an integer vector parallel
  to \code{dim(x)} that reports how the dimensions of \code{x} are mapped
  to the dimensions of its seed. If there is no \emph{net dimension
  rearrangement}, that is, if \code{dim(x)} is parallel to
  \code{dim(seed(x))}, the list has no \code{dimmap} attribute.

  If \code{as.DelayedOp} is \code{TRUE}, \code{summarizeMappingToSeed(x)}
  returns a linear tree with 2 \link{DelayedOp} nodes and a leaf node. The
  leaf node is the seed. Walking the tree from the seed, the 2 \link{DelayedOp}
  nodes are a \link{DelayedSubset} node followed by a \link{DelayedAperm}
  node. The former represents the \emph{net subsetting} and the latter the
  \emph{net dimension rearrangement}. Each of them can be a no-op.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedOp} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## showtree()
## ---------------------------------------------------------------------
m1 <- matrix(runif(150), nrow=15, ncol=10)
M1 <- DelayedArray(m1)

## By default, the tree of delayed operations carried by a DelayedArray
## object gets simplified each time a delayed operation is added to it.
## This can be disabled via a global option:
options(DelayedArray.simplify=FALSE)
M2 <- log(t(M1[5:1, c(TRUE, FALSE)] + 10))[-1, ]
showtree(M2)

## Note that as part of the simplification process, some operations
## can be reordered:
options(DelayedArray.simplify=TRUE)
M2 <- log(t(M1[5:1, c(TRUE, FALSE)] + 10))[-1, ]
showtree(M2)

## In the above example, the tree is linear i.e. all the operations
## are represented by unary nodes. The simplest way to know if a
## tree is linear is by counting its leaves with nseed():
nseed(M2)  # only 1 leaf means the tree is linear

options(DelayedArray.simplify=FALSE)

dimnames(M1) <- list(letters[1:15], LETTERS[1:10])
showtree(M1)

m2 <- matrix(1:20, nrow=10)
Y <- cbind(t(M1[ , 10:1]), DelayedArray(m2), M1[6:15, "A", drop=FALSE])
showtree(Y)
showtree(Y, show.node.dim=FALSE)
nseed(Y)  # the tree is not linear

Z <- t(Y[10:1, ])[1:15, ] + 0.4 * M1
showtree(Z)
nseed(Z)

Z@seed@seeds
Z@seed@seeds[[2]]@seed            # reaching to M1
Z@seed@seeds[[1]]@seed@seed@seed  # reaching to Y

## ---------------------------------------------------------------------
## summarizeMappingToSeed()
## ---------------------------------------------------------------------
a <- array(1:120, c(4, 5, 2))
M <- aperm(DelayedArray(a)[ , -1, ] / 100)[ , , 3] + 99:98
M
showtree(M)

summarizeMappingToSeed(M)  # 1st dimension was dropped, 2nd and 3rd
                           # dimension were permuted (transposition)

op2 <- summarizeMappingToSeed(M, as.DelayedOp=TRUE)
op2                # 2 nested delayed operations
op1 <- op2@seed
class(op1)         # DelayedSubset
class(op2)         # DelayedAperm
op1@index
op2@perm

DelayedArray(op2)  # same as M from a [, drop(), and aperm() point of
                   # view but the individual array elements are now reset
                   # to their original values i.e. to the values they
                   # have in the seed

## A simple function that returns TRUE if a DelayedArray object carries
## no "net subsetting" and no "net dimension rearrangement":
is_aligned_with_seed <- function(x)
{
    if (nseed(x) != 1L)
        return(FALSE)
    op2 <- summarizeMappingToSeed(x, as.DelayedOp=TRUE)
    op1 <- op2@seed
    isNoOp(op1) && isNoOp(op2)
}

M <- DelayedArray(a[ , , 1])
is_aligned_with_seed(log(M + 11:14) > 3)            # TRUE
is_aligned_with_seed(M[4:1, ])                      # FALSE
is_aligned_with_seed(M[4:1, ][4:1, ])               # TRUE
is_aligned_with_seed(t(M))                          # FALSE
is_aligned_with_seed(t(t(M)))                       # TRUE
is_aligned_with_seed(t(0.5 * t(M[4:1, ])[ , 4:1]))  # TRUE

options(DelayedArray.simplify=TRUE)
}

\keyword{methods}
