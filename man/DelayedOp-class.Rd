\name{DelayedOp-class}

% DelayedOp objects
\alias{class:DelayedOp}
\alias{DelayedOp-class}
\alias{DelayedOp}

\alias{path,DelayedOp-method}
\alias{path<-,DelayedOp-method}

% DelayedSubset objects
\alias{dim,DelayedSubset-method}
\alias{dimnames,DelayedSubset-method}
\alias{extract_array,DelayedSubset-method}

% DelayedDimnames objects
\alias{dim,DelayedDimnames-method}
\alias{dimnames,DelayedDimnames-method}
\alias{extract_array,DelayedDimnames-method}

% DelayedUnaryIsoOp objects
\alias{dim,DelayedUnaryIsoOp-method}
\alias{dimnames,DelayedUnaryIsoOp-method}
\alias{extract_array,DelayedUnaryIsoOp-method}

% DelayedAperm objects
\alias{dim,DelayedAperm-method}
\alias{dimnames,DelayedAperm-method}
\alias{extract_array,DelayedAperm-method}

% DelayedVariadicIsoOp objects
\alias{dim,DelayedVariadicIsoOp-method}
\alias{dimnames,DelayedVariadicIsoOp-method}
\alias{extract_array,DelayedVariadicIsoOp-method}

% DelayedAbind objects
\alias{dim,DelayedAbind-method}
\alias{dimnames,DelayedAbind-method}
\alias{extract_array,DelayedAbind-method}

\title{DelayedOp objects}

\description{
  In a \link{DelayedArray} object the delayed operations are stored as a tree
  of DelayedOp objects. Each node in the tree is represented by a DelayedOp
  object.

  DelayedOp objects are used inside \link{DelayedArray} objects and are not
  intended to be manipulated directly by the end user.
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
  Unary nodes (i.e. nodes with an outdegree of 1) must also support the
  \code{seed()} and \code{path()} getters and setters.
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

\keyword{internal}
