\name{DelayedOp-class}

% DelayedOp objects
\alias{class:DelayedOp}
\alias{DelayedOp-class}
\alias{DelayedOp}

\alias{showtree}
\alias{showtree,ANY-method}

% DelayedSubset objects
\alias{summary.DelayedSubset}
\alias{summary,DelayedSubset-method}
\alias{dim,DelayedSubset-method}
\alias{dimnames,DelayedSubset-method}
\alias{extract_array,DelayedSubset-method}

% DelayedDimnames objects
\alias{summary.DelayedDimnames}
\alias{summary,DelayedDimnames-method}
\alias{dim,DelayedDimnames-method}
\alias{dimnames,DelayedDimnames-method}
\alias{extract_array,DelayedDimnames-method}

% DelayedUnaryIsoOp objects
\alias{summary.DelayedUnaryIsoOp}
\alias{summary,DelayedUnaryIsoOp-method}
\alias{dim,DelayedUnaryIsoOp-method}
\alias{dimnames,DelayedUnaryIsoOp-method}
\alias{extract_array,DelayedUnaryIsoOp-method}

% DelayedAperm objects
\alias{summary.DelayedAperm}
\alias{summary,DelayedAperm-method}
\alias{dim,DelayedAperm-method}
\alias{dimnames,DelayedAperm-method}
\alias{extract_array,DelayedAperm-method}

% DelayedVariadicIsoOp objects
\alias{summary.DelayedVariadicIsoOp}
\alias{summary,DelayedVariadicIsoOp-method}
\alias{dim,DelayedVariadicIsoOp-method}
\alias{dimnames,DelayedVariadicIsoOp-method}
\alias{extract_array,DelayedVariadicIsoOp-method}

% DelayedAbind objects
\alias{summary.DelayedAbind}
\alias{summary,DelayedAbind-method}
\alias{dim,DelayedAbind-method}
\alias{dimnames,DelayedAbind-method}
\alias{extract_array,DelayedAbind-method}

% updateObject()
\alias{updateObject,DelayedOp-method}
\alias{updateObject,SeedDimPicker-method}
\alias{updateObject,ConformableSeedCombiner-method}
\alias{updateObject,SeedBinder-method}

\title{DelayedOp objects}

\description{
  In a \link{DelayedArray} object the delayed operations are stored as a tree
  of DelayedOp objects. Each node in the tree is represented by a DelayedOp
  object.

  DelayedOp objects are used inside \link{DelayedArray} objects and are not
  intended to be manipulated directly by the end user.

  \code{showtree} can be used to visualize the tree of delayed operations
  in a \link{DelayedArray} object.
}

\usage{
showtree(x)
}

\arguments{
  \item{x}{
    Typically a \link{DelayedArray} object. Can also be a DelayedOp object
    or a list where each element is a \link{DelayedArray} or DelayedOp object.
  }
}

\details{
  6 types of nodes are currently supported. Each type is a DelayedOp subclass:

  \preformatted{  Node type      Outdegree  Operation
  ---------------------------------------------------------------------
  DelayedSubset          1  Multi-dimensional single bracket subsetting
  DelayedDimnames        1  Set dimnames
  DelayedUnaryIsoOp      1  Unary op that preserves the geometry
  DelayedAperm           1  Extended aperm() (can drop dimensions)
  DelayedVariadicIsoOp   N  N-ary op that preserves the geometry
  DelayedAbind           N  abind()
  }

  All the nodes are array-like objects that must satisfy the \emph{seed
  contract} i.e. they must support \code{dim()}, \code{dimnames()}, and
  \code{extract_array()}. See \code{?\link{extract_array}} for more
  information about the \emph{seed contract}.
}

\note{
  The DelayedOp virtual class and its 6 concrete subclasses are for
  internal use only and are not exported.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \code{\link{extract_array}}.
  }
}

\examples{
m1 <- matrix(runif(150), nrow=15, ncol=10)
X <- DelayedArray(m1)
showtree(X)

dimnames(X) <- list(letters[1:15], LETTERS[1:10])
showtree(X)

m2 <- matrix(1:20, nrow=10)
Y <- cbind(t(X[ , 10:1]), DelayedArray(m2), X[6:15, "A"])
showtree(Y)

Z <- t(Y[10:1, ])[1:15, ] + 0.4 * X
showtree(Z)
showtree(Z@seed@seeds)
showtree(Z@seed@seeds[[2]]@seed)  # reaching to X
showtree(Z@seed@seeds[[1]]@seed@seed@seed)  # reaching to Y
}

\keyword{internal}
