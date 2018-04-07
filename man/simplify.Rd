\name{simplify}

% showtree()
\alias{showtree}
\alias{showtree,ANY-method}
\alias{show,DelayedOp-method}

% simplify()
\alias{simplify}
\alias{simplify,DelayedSubset-method}
\alias{simplify,DelayedAperm-method}
\alias{simplify,DelayedUnaryIsoOp-method}
\alias{simplify,DelayedDimnames-method}
\alias{simplify,DelayedArray-method}

\title{Simplify a tree of delayed operations}

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
options(DelayedArray.simplify=FALSE)

m1 <- matrix(runif(150), nrow=15, ncol=10)
X <- DelayedArray(m1)
showtree(X)
nseed(X)  # the tree is linear

dimnames(X) <- list(letters[1:15], LETTERS[1:10])
showtree(X)

m2 <- matrix(1:20, nrow=10)
Y <- cbind(t(X[ , 10:1]), DelayedArray(m2), X[6:15, "A", drop=FALSE])
showtree(Y)
showtree(Y, show.node.dim=FALSE)
nseed(Y)  # the tree is not linear

Z <- t(Y[10:1, ])[1:15, ] + 0.4 * X
showtree(Z)
nseed(Z)

Z@seed@seeds
Z@seed@seeds[[2]]@seed            # reaching to X
Z@seed@seeds[[1]]@seed@seed@seed  # reaching to Y

options(DelayedArray.simplify=TRUE)
}

\keyword{methods}
