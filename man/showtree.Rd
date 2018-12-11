\name{showtree}

\alias{showtree}
\alias{show,DelayedOp-method}

\alias{nseed}
\alias{nseed,ANY-method}

\alias{seed}
\alias{seed,DelayedOp-method}
\alias{seed<-}
\alias{seed<-,DelayedOp-method}

\alias{path}
\alias{path,DelayedOp-method}
\alias{path<-,DelayedOp-method}

\alias{seedApply}

\title{Visualize and access the leaves of a tree of delayed operations}

\description{
  \code{showtree} can be used to visualize the tree of delayed operations
  carried by a \link{DelayedArray} object.

  Use \code{nseed}, \code{seed}, or \code{path} to access the number of seeds,
  the seed, or the seed path of a \link{DelayedArray} object, respectively.

  Use \code{seedApply} to apply a function to the seeds of a
  \link{DelayedArray} object.
}

\usage{
showtree(x, show.node.dim=TRUE)

nseed(x)            # seed counter
seed(x)             # seed getter and setter
path(object, ...)   # path getter and setter

seedApply(x, FUN, ...)
}

\arguments{
  \item{x, object}{
    Typically a \link{DelayedArray} object but can also be a \link{DelayedOp}
    object or a list where each element is a \link{DelayedArray} or
    \link{DelayedOp} object.
  }
  \item{show.node.dim}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} (the default), the nodes
    dimensions and data type are displayed.
  }
  \item{FUN}{
    The function to be applied to each leaf in \code{x}.
  }
  \item{...}{
    Optional arguments to \code{FUN} for \code{seedApply()}.

    Additional arguments passed to methods for \code{path()}.
  }
}

\value{
  The number of seeds contained in \code{x} for \code{nseed}.

  The seed contained in \code{x} for \code{seed}.

  The path of the seed contained in \code{object} for \code{path}.

  A list of length \code{nseed(x)} for \code{seedApply}.
}

\seealso{
  \itemize{
    \item \code{\link{simplify}} to simplify the tree of delayed operations
          carried by a \link{DelayedArray} object.

    \item \link{DelayedOp} objects.

    \item \link{DelayedArray} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## showtree(), nseed(), and seed()
## ---------------------------------------------------------------------
m1 <- matrix(runif(150), nrow=15, ncol=10)
M1 <- DelayedArray(m1)
showtree(M1)
seed(M1)

M2 <- log(t(M1[5:1, c(TRUE, FALSE)] + 10))[-1, ]
showtree(M2)

## In the above example, the tree is linear i.e. all the operations
## are represented by unary nodes. The simplest way to know if a
## tree is linear is by counting its leaves with nseed():
nseed(M2)  # only 1 leaf means the tree is linear
seed(M2)

dimnames(M1) <- list(letters[1:15], LETTERS[1:10])
showtree(M1)

m2 <- matrix(1:20, nrow=10)
Y <- cbind(t(M1[ , 10:1]), DelayedArray(m2), M1[6:15, "A", drop=FALSE])
showtree(Y)
showtree(Y, show.node.dim=FALSE)
nseed(Y)  # the tree is not linear

Z <- t(Y[10:1, ])[1:15, ] + 0.4 * M1
showtree(Z)
nseed(Z)  # the tree is not linear

## ---------------------------------------------------------------------
## seedApply()
## ---------------------------------------------------------------------
seedApply(Y, class)
seedApply(Y, dim)
}

\keyword{methods}
