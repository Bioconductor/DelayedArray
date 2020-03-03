\name{bind-arrays}

\alias{bind-arrays}
\alias{bind arrays}

\alias{arbind}
\alias{acbind}
\alias{arbind,array-method}
\alias{acbind,array-method}


\title{Bind arrays along their rows or columns}

\description{
  Bind array-like objects with an arbitrary number of dimensions along their
  rows (\code{arbind}) or columns (\code{acbind}).
}

\usage{
arbind(...)
acbind(...)
}

\arguments{
  \item{...}{
    The array-like objects to bind.
  }
}

\value{
  An array-like object, typically of the same class as the input objects if
  they all have the same class.
}

\seealso{
  \itemize{
    \item \code{\link{DelayedArray}} in this package for
          arbind/acbind'ing DelayedArray objects.

    \item \code{\link[base]{rbind}} and \code{\link[base]{cbind}} in the
          \pkg{base} package for the corresponding operations on matrix-like
          objects.

    \item The \pkg{abind} package from CRAN.
  }
}

\examples{
a1 <- array(1:60, c(3, 5, 4),
            dimnames=list(NULL, paste0("M1y", 1:5), NULL))
a2 <- array(101:240, c(7, 5, 4),
            dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL))
a3 <- array(10001:10100, c(5, 5, 4),
            dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))

arbind(a1, a2, a3)
}
\keyword{methods}
