\name{DelayedOp-class}

% DelayedOp objects
\alias{class:DelayedOp}
\alias{DelayedOp-class}
\alias{DelayedOp}

\alias{isNoOp}

\alias{summary.DelayedOp}
\alias{summary,DelayedOp-method}

% DelayedUnaryOp objects
\alias{class:DelayedUnaryOp}
\alias{DelayedUnaryOp-class}
\alias{DelayedUnaryOp}

\alias{isSparse,DelayedUnaryOp-method}

% DelayedSubset objects
\alias{class:DelayedSubset}
\alias{DelayedSubset-class}
\alias{DelayedSubset}

\alias{isNoOp,DelayedSubset-method}
\alias{summary.DelayedSubset}
\alias{summary,DelayedSubset-method}
\alias{dim,DelayedSubset-method}
\alias{dimnames,DelayedSubset-method}
\alias{extract_array,DelayedSubset-method}
\alias{isSparse,DelayedSubset-method}

% DelayedAperm objects
\alias{class:DelayedAperm}
\alias{DelayedAperm-class}
\alias{DelayedAperm}

\alias{isNoOp,DelayedAperm-method}
\alias{summary.DelayedAperm}
\alias{summary,DelayedAperm-method}
\alias{dim,DelayedAperm-method}
\alias{dimnames,DelayedAperm-method}
\alias{extract_array,DelayedAperm-method}

% DelayedUnaryIsoOp objects
\alias{class:DelayedUnaryIsoOp}
\alias{DelayedUnaryIsoOp-class}
\alias{DelayedUnaryIsoOp}

\alias{dim,DelayedUnaryIsoOp-method}
\alias{dimnames,DelayedUnaryIsoOp-method}
\alias{extract_array,DelayedUnaryIsoOp-method}

% DelayedUnaryIsoOpStack objects
\alias{class:DelayedUnaryIsoOpStack}
\alias{DelayedUnaryIsoOpStack-class}
\alias{DelayedUnaryIsoOpStack}

\alias{summary.DelayedUnaryIsoOpStack}
\alias{summary,DelayedUnaryIsoOpStack-method}
\alias{extract_array,DelayedUnaryIsoOpStack-method}
\alias{isSparse,DelayedUnaryIsoOpStack-method}

% DelayedUnaryIsoOpWithArgs objects
\alias{class:DelayedUnaryIsoOpWithArgs}
\alias{DelayedUnaryIsoOpWithArgs-class}
\alias{DelayedUnaryIsoOpWithArgs}

\alias{summary.DelayedUnaryIsoOpWithArgs}
\alias{summary,DelayedUnaryIsoOpWithArgs-method}
\alias{extract_array,DelayedUnaryIsoOpWithArgs-method}
\alias{isSparse,DelayedUnaryIsoOpWithArgs-method}

% DelayedDimnames objects
\alias{class:DelayedDimnames}
\alias{DelayedDimnames-class}
\alias{DelayedDimnames}

\alias{isNoOp,DelayedDimnames-method}
\alias{summary.DelayedDimnames}
\alias{summary,DelayedDimnames-method}
\alias{dimnames,DelayedDimnames-method}

% DelayedNaryOp objects
\alias{class:DelayedNaryOp}
\alias{DelayedNaryOp-class}
\alias{DelayedNaryOp}

% DelayedNaryIsoOp objects
\alias{summary.DelayedNaryIsoOp}
\alias{summary,DelayedNaryIsoOp-method}
\alias{dim,DelayedNaryIsoOp-method}
\alias{dimnames,DelayedNaryIsoOp-method}
\alias{extract_array,DelayedNaryIsoOp-method}
\alias{isSparse,DelayedNaryIsoOp-method}

% DelayedAbind objects
\alias{summary.DelayedAbind}
\alias{summary,DelayedAbind-method}
\alias{dim,DelayedAbind-method}
\alias{dimnames,DelayedAbind-method}
\alias{extract_array,DelayedAbind-method}
\alias{isSparse,DelayedAbind-method}

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

  \code{\link{showtree}} and \code{\link{simplify}} can be used to visualize
  and simplify this tree.
}

\usage{
isNoOp(x)
}

\arguments{
  \item{x}{
    A DelayedSubset, DelayedAperm, or DelayedDimnames object.
  }
}

\details{
  7 types of nodes are currently supported. Each type is a DelayedOp subclass:

  \preformatted{  Node type                        Represented operation
  ------------------------------------------------------------------
  DelayedOp (VIRTUAL)
  ------------------------------------------------------------------
    DelayedUnaryOp (VIRTUAL)
    - DelayedSubset                Multi-dimensional single bracket
                                   subsetting.
    - DelayedAperm                 Extended aperm() (can drop
                                   dimensions).
    - DelayedUnaryIsoOp (VIRTUAL)  Unary op that preserves the
                                   geometry.
      - DelayedUnaryIsoOpStack     Simple ops stacked together.
      - DelayedUnaryIsoOpWithArgs  One op with vector-like arguments
                                   along the dimensions of the input.
      - DelayedDimnames            Set/replace the dimnames.
  ------------------------------------------------------------------
    DelayedNaryOp (VIRTUAL)
    - DelayedNaryIsoOp             N-ary op that preserves the
                                   geometry.
    - DelayedAbind                 abind()
  ------------------------------------------------------------------
  }

  All the nodes are array-like objects that must comply with the \emph{seed
  contract} i.e. they must support \code{dim()}, \code{dimnames()}, and
  \code{extract_array()}. See \code{?\link{extract_array}} for more
  information about the \emph{seed contract}.

  \code{isNoOp()} can be called on a DelayedSubset, DelayedAperm, or
  DelayedDimnames object and will return \code{TRUE} if the object
  represents a no-op.
}

\note{
  The DelayedOp virtual class and its 7 concrete subclasses are for
  internal use only and are not exported.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \code{\link{showtree}} to visualize, simplify, and inspect the
          tree of delayed operations in a \link{DelayedArray} object.

    \item \code{\link{extract_array}}.
  }
}

\keyword{methods}
