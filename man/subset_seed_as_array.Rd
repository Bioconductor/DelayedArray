\name{subset_seed_as_array}

\alias{subset_seed_as_array}
\alias{subset_seed_as_array,ANY-method}
\alias{subset_seed_as_array,array-method}
\alias{subset_seed_as_array,data.frame-method}
\alias{subset_seed_as_array,DataFrame-method}

\alias{dim,SeedBinder-method}
\alias{dimnames,SeedBinder-method}
\alias{subset_seed_as_array,SeedBinder-method}

\title{subset_seed_as_array}

\description{
  \code{subset_seed_as_array} is an internal generic function not
  aimed to be used directly by the user.
  It has methods defined for array, data.frame, DataFrame objects
  and other array-like objects.

  The \code{DelayedArray()} constructor function will accept any
  seed that supports \code{dim()}, \code{dimnames()}, and
  \code{subset_seed_as_array()}.
}

\usage{
subset_seed_as_array(seed, index)
}

\arguments{
  \item{seed}{
    An array-like object.
  }
  \item{index}{
    An unnamed list of subscripts as positive integer vectors, one vector
    per seed dimension. \emph{Missing} list elements are allowed and must
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
