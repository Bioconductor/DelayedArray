\name{extract_array}

\alias{extract_array}
\alias{extract_array,ANY-method}
\alias{extract_array,array-method}
\alias{extract_array,data.frame-method}
\alias{extract_array,DataFrame-method}

\alias{type}
\alias{type,array-method}
\alias{type,ANY-method}

\title{extract_array}

\description{
  \code{extract_array} is an internal generic function not aimed to be used
  directly by the user. It has methods defined for array, data.frame,
  \link[S4Vectors]{DataFrame} objects and other array-like objects.

  The \code{DelayedArray()} constructor function will accept any seed that
  satisfies the \emph{seed contract} i.e. that supports \code{dim()},
  \code{dimnames()}, and \code{extract_array()}.
}

\usage{
extract_array(x, index)

type(x)
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{index}{
    An unnamed list of subscripts as positive integer vectors, one vector
    per dimension in \code{x}. Empty and missing subscripts (represented
    by \code{integer(0)} and \code{NULL} list elements, respectively) are
    allowed. The subscripts can contain duplicated indices. They cannot
    contain NAs or non-positive values.
  }
}

\details{
  \code{extract_array} methods need to support empty and missing subscripts
  e.g. \code{extract_array(x, list(NULL, integer(0)))} must return an M x 0
  matrix and \code{extract_array(x, list(integer(0), integer(0)))} a 0 x 0
  matrix.

  Also subscripts are allowed to contain duplicated indices so things like
  \code{extract_array(seed, list(c(1:3, 3:1), 2L))} need to be supported.

  \code{type(x)} returns the type of the elements in \code{x}. It's equivalent
  to \code{typeof(x)} or \code{storage.mode(x)} on an ordinary array.
  It works out-of-the-box on any array-like object \code{x} for which
  \code{extract_array(x)} works.
}

\value{
  \code{extract_array} must return an \emph{ordinary} array of the
  appropriate type (i.e. integer, double, etc...). For example, if
  \code{x} is an object representing an M x N matrix of complex numbers,
  \code{extract_array(x, list(NULL, 2L))} must return its 2nd column as
  an \emph{ordinary} M x 1 matrix of type complex.

  \code{type} returns the type of the array elements.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link[base]{array} and \link[base]{data.frame} objects in base R.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.
  }
}

\keyword{internal}
