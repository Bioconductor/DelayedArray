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

% DelayedSubset objects
\alias{class:DelayedSubset}
\alias{DelayedSubset-class}
\alias{DelayedSubset}

\alias{is_noop,DelayedSubset-method}
\alias{summary.DelayedSubset}
\alias{summary,DelayedSubset-method}

\alias{dim,DelayedSubset-method}
\alias{dimnames,DelayedSubset-method}
\alias{extract_array,DelayedSubset-method}

\alias{is_sparse,DelayedSubset-method}
\alias{extract_sparse_array,DelayedSubset-method}

% DelayedAperm objects
\alias{class:DelayedAperm}
\alias{DelayedAperm-class}
\alias{DelayedAperm}

\alias{is_noop,DelayedAperm-method}
\alias{summary.DelayedAperm}
\alias{summary,DelayedAperm-method}

\alias{dim,DelayedAperm-method}
\alias{dimnames,DelayedAperm-method}
\alias{extract_array,DelayedAperm-method}

\alias{is_sparse,DelayedAperm-method}
\alias{extract_sparse_array,DelayedAperm-method}

% DelayedUnaryIsoOp objects
\alias{class:DelayedUnaryIsoOp}
\alias{DelayedUnaryIsoOp-class}
\alias{DelayedUnaryIsoOp}

\alias{dim,DelayedUnaryIsoOp-method}
\alias{dimnames,DelayedUnaryIsoOp-method}
\alias{extract_array,DelayedUnaryIsoOp-method}

\alias{is_sparse,DelayedUnaryIsoOp-method}
\alias{extract_sparse_array,DelayedUnaryIsoOp-method}

% DelayedUnaryIsoOpStack objects
\alias{class:DelayedUnaryIsoOpStack}
\alias{DelayedUnaryIsoOpStack-class}
\alias{DelayedUnaryIsoOpStack}

\alias{summary.DelayedUnaryIsoOpStack}
\alias{summary,DelayedUnaryIsoOpStack-method}

\alias{extract_array,DelayedUnaryIsoOpStack-method}

\alias{is_sparse,DelayedUnaryIsoOpStack-method}
\alias{extract_sparse_array,DelayedUnaryIsoOpStack-method}

% DelayedUnaryIsoOpWithArgs objects
\alias{class:DelayedUnaryIsoOpWithArgs}
\alias{DelayedUnaryIsoOpWithArgs-class}
\alias{DelayedUnaryIsoOpWithArgs}

\alias{summary.DelayedUnaryIsoOpWithArgs}
\alias{summary,DelayedUnaryIsoOpWithArgs-method}

\alias{extract_array,DelayedUnaryIsoOpWithArgs-method}

\alias{is_sparse,DelayedUnaryIsoOpWithArgs-method}
\alias{extract_sparse_array,DelayedUnaryIsoOpWithArgs-method}

% DelayedSubassign objects
\alias{class:DelayedSubassign}
\alias{DelayedSubassign-class}
\alias{DelayedSubassign}

\alias{is_noop,DelayedSubassign-method}
\alias{summary.DelayedSubassign}
\alias{summary,DelayedSubassign-method}

\alias{extract_array,DelayedSubassign-method}

\alias{is_sparse,DelayedSubassign-method}
\alias{extract_sparse_array,DelayedSubassign-method}

% DelayedDimnames objects
\alias{class:DelayedDimnames}
\alias{DelayedDimnames-class}
\alias{DelayedDimnames}

\alias{is_noop,DelayedDimnames-method}
\alias{summary.DelayedDimnames}
\alias{summary,DelayedDimnames-method}

\alias{dimnames,DelayedDimnames-method}

% DelayedNaryOp objects
\alias{class:DelayedNaryOp}
\alias{DelayedNaryOp-class}
\alias{DelayedNaryOp}

% DelayedNaryIsoOp objects
\alias{class:DelayedNaryIsoOp}
\alias{DelayedNaryIsoOp-class}
\alias{DelayedNaryIsoOp}

\alias{summary.DelayedNaryIsoOp}
\alias{summary,DelayedNaryIsoOp-method}

\alias{dim,DelayedNaryIsoOp-method}
\alias{dimnames,DelayedNaryIsoOp-method}
\alias{extract_array,DelayedNaryIsoOp-method}

\alias{is_sparse,DelayedNaryIsoOp-method}
\alias{extract_sparse_array,DelayedNaryIsoOp-method}

% DelayedAbind objects
\alias{class:DelayedAbind}
\alias{DelayedAbind-class}
\alias{DelayedAbind}

\alias{summary.DelayedAbind}
\alias{summary,DelayedAbind-method}

\alias{dim,DelayedAbind-method}
\alias{dimnames,DelayedAbind-method}
\alias{extract_array,DelayedAbind-method}

\alias{is_sparse,DelayedAbind-method}
\alias{extract_sparse_array,DelayedAbind-method}

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
is_noop(x)
}

\arguments{
  \item{x}{
    A DelayedSubset, DelayedAperm, or DelayedDimnames object.
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
      - DelayedDimnames            Set/replace the dimnames.
  ------------------------------------------------------------------
  * DelayedNaryOp (VIRTUAL)
    o DelayedNaryIsoOp             N-ary op that preserves the
                                   geometry.
    o DelayedAbind                 abind()
  ------------------------------------------------------------------
  }

  All the nodes are array-like objects that must comply with the \emph{seed
  contract} i.e. they must support \code{dim()}, \code{dimnames()}, and
  \code{extract_array()}. See \code{?\link{extract_array}} for more
  information about the \emph{seed contract}.

  \code{is_noop()} can only be called on a DelayedSubset, DelayedAperm,
  or DelayedDimnames object at the moment, and will return \code{TRUE}
  if the object represents a no-op.
}

\note{
  The DelayedOp virtual class and its 8 concrete subclasses are for
  internal use only and never exposed to the end user.
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
