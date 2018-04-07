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

\title{Visualize and simplify a tree of delayed operations}

\description{
  \code{showtree} can be used to visualize the tree of delayed operations
  in a \link{DelayedArray} object.

  \code{simplify} can be used to simplify it.
}

\usage{
showtree(x, show.node.dim=TRUE)

simplify(x)
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
    dimensions and data type is displayed.
  }
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedOp} objects.
  }
}

\examples{
m1 <- matrix(runif(150), nrow=15, ncol=10)
M1 <- DelayedArray(m1)

## By default, the tree of delayed operations in a DelayedArray object
## gets simplified. This can be disabled via a global option:
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

options(DelayedArray.simplify=TRUE)
}

\keyword{methods}
