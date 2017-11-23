\name{extract_array}

\alias{extract_array}
\alias{extract_array,ANY-method}
\alias{extract_array,array-method}
\alias{extract_array,data.frame-method}
\alias{extract_array,DataFrame-method}

\alias{dim,ConformableSeedCombiner-method}
\alias{dimnames,ConformableSeedCombiner-method}
\alias{extract_array,ConformableSeedCombiner-method}

\alias{dim,SeedBinder-method}
\alias{dimnames,SeedBinder-method}
\alias{extract_array,SeedBinder-method}

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
}

\arguments{
  \item{x}{
    An array-like object.
  }
  \item{index}{
    An unnamed list of subscripts as positive integer vectors, one vector
    per dimension in \code{x}. \emph{Missing} list elements are allowed and must
    be represented by \code{NULL}s.
  }
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link[base]{array} and \link[base]{data.frame} objects in base R.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.
  }
}

\keyword{internal}
