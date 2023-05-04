\name{DelayedOp-class}

% DelayedOp objects
\alias{class:DelayedOp}
\alias{DelayedOp-class}
\alias{DelayedOp}

\alias{is_noop}
\alias{summary.DelayedOp}
\alias{summary,DelayedOp-method}

% DelayedUnaryOp objects
\alias{class:DelayedUnaryOp}
\alias{DelayedUnaryOp-class}
\alias{DelayedUnaryOp}

% DelayedUnaryIsoOp objects
\alias{class:DelayedUnaryIsoOp}
\alias{DelayedUnaryIsoOp-class}
\alias{DelayedUnaryIsoOp}

\alias{dim,DelayedUnaryIsoOp-method}
\alias{dimnames,DelayedUnaryIsoOp-method}
\alias{extract_array,DelayedUnaryIsoOp-method}

\alias{is_sparse,DelayedUnaryIsoOp-method}
\alias{OLD_extract_sparse_array,DelayedUnaryIsoOp-method}

% DelayedNaryOp objects
\alias{class:DelayedNaryOp}
\alias{DelayedNaryOp-class}
\alias{DelayedNaryOp}

% updateObject()
\alias{updateObject,DelayedOp-method}

\title{DelayedOp objects}

\description{
  NOTE: This man page is about \link{DelayedArray} internals and is provided
  for developers and advanced users only.

  In a \link{DelayedArray} object, the delayed operations are stored as a
  tree where the leaves are operands and the nodes are the operations.
  Each node in the tree is a DelayedOp derivative representing a particular
  delayed operation.

  DelayedOp is a virtual class with 8 concrete subclasses. Each subclass
  provides a formal representation for a particular kind of delayed
  operation.
}

\usage{
is_noop(x)
}

\arguments{
  \item{x}{
    A DelayedSubset, DelayedAperm, or DelayedSetDimnames object.
  }
}

\details{
  8 types of nodes are currently supported. Each type is a DelayedOp subclass:

  \preformatted{  Node type                        Represented operation
  ------------------------------------------------------------------
  DelayedOp (VIRTUAL)
  ------------------------------------------------------------------
  * DelayedUnaryOp (VIRTUAL)
    o DelayedSubset                Multi-dimensional single bracket
                                   subsetting.
    o DelayedAperm                 Extended aperm() (can drop and/or
                                   add ineffective dimensions).
    o DelayedUnaryIsoOp (VIRTUAL)  Unary op that preserves the
                                   geometry.
      - DelayedUnaryIsoOpStack     Simple ops stacked together.
      - DelayedUnaryIsoOpWithArgs  One op with vector-like arguments
                                   along the dimensions of the input.
      - DelayedSubassign           Multi-dimensional single bracket
                                   subassignment.
      - DelayedSetDimnames         Set/replace the dimnames.
  ------------------------------------------------------------------
  * DelayedNaryOp (VIRTUAL)
    o DelayedNaryIsoOp             N-ary op that preserves the
                                   geometry.
    o DelayedAbind                 abind()
  ------------------------------------------------------------------
  }

  All the nodes are array-like objects that must comply with the \emph{seed
  contract} i.e. they must support \code{dim()}, \code{dimnames()}, and
  \code{extract_array()}. See \code{?\link[S4Arrays]{extract_array}}
  in the \pkg{S4Arrays} package for more information about the
  \emph{seed contract}.

  \code{is_noop()} can only be called on a DelayedSubset, DelayedAperm,
  or DelayedSetDimnames object at the moment, and will return \code{TRUE}
  if the object represents a no-op.
}

\note{
  The DelayedOp virtual class and its 8 concrete subclasses are used inside
  a \link{DelayedArray} object to represent delayed operations carried by
  the object. They're never exposed to the end user and are not intended to
  be manipulated directly.
}

\seealso{
  \itemize{
    \item DelayedOp concrete subclasses: \link{DelayedSubset},
          \link{DelayedAperm}, \link{DelayedUnaryIsoOpStack},
          \link{DelayedUnaryIsoOpWithArgs}, \link{DelayedSubassign},
          \link{DelayedSetDimnames}, \link{DelayedNaryIsoOp},
          and \link{DelayedAbind}.

    \item \link{DelayedArray} objects.

    \item \code{\link{showtree}} to visualize the nodes and access the
          leaves in the tree of delayed operations carried by a
          \link{DelayedArray} object.

    \item \code{\link{simplify}} to simplify the tree of delayed
          operations carried by a \link{DelayedArray} object.

    \item \code{\link[S4Arrays]{extract_array}} in the \pkg{S4Arrays} package.
  }
}

\keyword{methods}
