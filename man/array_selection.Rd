\name{array selection}

\alias{array selection}
\alias{array_selection}
\alias{array selections}
\alias{array_selections}
\alias{L-index}
\alias{Lindex}
\alias{M-index}
\alias{Mindex}
\alias{N-index}
\alias{Nindex}

\alias{Lindex2Mindex}
\alias{Mindex2Lindex}
\alias{linearInd}

\title{Manipulation of array selections}

\description{
  NOTE: The tools documented in this man page are primarily intended
  for developers or advanced users curious about the internals of the
  \pkg{DelayedArray} package. End users typically don't need them for
  their regular use of \link{DelayedArray} objects.

  An \emph{array selection} is just an index into an array-like object that
  contains the information of which array elements are selected. This index
  can take various forms but 3 special forms are particularly useful and
  extensively used thoughout the \pkg{DelayedArray} framework:
  \emph{linear index} (also referred to as \emph{L-index} or \emph{Lindex}),
  \emph{matrix index} (also referred to as \emph{M-index} or \emph{Mindex}),
  \emph{N-dimensional index} (also referred to as \emph{N-index} or
  \emph{Nindex}). See Details section below for more information.

  Two utility functions are provided at the moment to convert back and forth
  between \emph{L-indices} and \emph{M-indices}. More will be added in the
  future to convert between other types of array indices.
}

\usage{
## Convert back and forth between L-indices and M-indices:
Lindex2Mindex(Lindex, dim, use.names=FALSE)
Mindex2Lindex(Mindex, dim, use.names=FALSE, as.integer=FALSE)
}

\arguments{
  \item{Lindex}{
    An \emph{L-index}. See Details section below.
  }
  \item{Mindex}{
    An \emph{M-index}. See Details section below.
    For convenience, \code{Mindex} can also be specified as an integer vector
    with one element per dimension in the underlying array, in which case it
    will be treated like a 1-row matrix.
  }
  \item{dim}{
    An integer vector containing the dimensions of the underlying array.

    Note that \code{dim} can also be an integer matrix, in which case
    it must have one row per element in \code{Lindex} (or per row in
    \code{Mindex}) and one column per dimension in the underlying array.
  }
  \item{use.names}{
    Should the names (or rownames) on the input be propagated to the output?
  }
  \item{as.integer}{
    Set to \code{TRUE} to force \code{Mindex2Lindex} to return the L-index
    as an integer vector. Dangerous!

    By default, i.e. when \code{as.integer=FALSE}, \code{Mindex2Lindex} will
    return the L-index either as an integer or numeric vector. It will choose
    the former only if it's safe, that is, only if all the values in the
    L-index "fit" in the integer type. More precisely:
    \itemize{
      \item If \code{dim} is not a matrix (i.e. is a vector) or if it's a
            matrix with a single row: \code{Mindex2Lindex} returns an integer
            or numeric vector depending on whether \code{prod(dim)} is <=
            \code{.Machine$integer.max} (2^31 - 1) or not.
      \item Otherwise \code{Mindex2Lindex} returns a numeric vector.
    }
    Note that with these rules, \code{Mindex2Lindex} can return a numeric
    vector even if an integer vector could have been used.

    Use \code{as.integer=TRUE} only in situations where you know that all
    the L-index values are going to "fit" in the integer type.
    \code{Mindex2Lindex} will return garbage if they don't.
  }
}

\details{
  The 3 special forms of array indices extensively used thoughout the
  \pkg{DelayedArray} framework:

  \enumerate{
    \item \emph{Linear index} (or \emph{L-index} or \emph{Lindex}):
          A numeric vector with no NAs where each value is >= 1 and <=
          the length of the array-like object.
          When using an L-index to subset an array-like object, the returned
          value is a vector-like object (i.e. no dimensions) of the same
          length as the L-index.

          Example:
          \preformatted{
            a <- array(101:124, 4:2)
            Lindex <- c(7, 2, 24, 2)
            a[Lindex]
          }

    \item \emph{Matrix index} (or \emph{M-index} or \emph{Mindex}):
          An integer matrix with one column per dimension in the array-like
          object and one row per array element in the selection. No NAs.
          The values in each column must be >= 1 and <= the extent of the
          array-like object along the corresponding dimension.
          When using an M-index to subset an array-like object, the returned
          value is a vector-like object (i.e. no dimensions) of length the
          number of rows in the M-index.

          Example:
          \preformatted{
            a <- array(101:124, 4:2)
            Mindex <- rbind(c(3, 2, 1),
                            c(2, 1, 1),
                            c(4, 3, 2),
                            c(2, 1, 1))
            a[Mindex]
          }

          Note that this is the type of index returned by
          \code{base::\link[base]{arrayInd}}.

    \item \emph{N-dimensional} (or \emph{N-index} or \emph{Nindex}):
          A list with one list element per dimension in the array-like object.
          Each list element must be a subscript describing the selection along
          the corresponding dimension of the array-like object.
          IMPORTANT: A \code{NULL} subscript is interpreted as a \emph{missing}
          subscript ("missing" like in \code{a[ , , 1:2]}), that is, as a
          subscript that runs along the full extend of the corresponding
          dimension of the array-like object. This means that before an
          N-index can be used in a call to \code{[}, \code{[<-}, \code{[[}
          or \code{[[<-}, the \code{NULL} list elements in it must be
          replaced with objects of class \code{"name"}.
          When using an N-index to subset an array-like object, the returned
          value is another array-like object of dimensions the lengths of the
          selections along each dimensions.

          Examples:
          \preformatted{
            a <- array(101:124, 4:2)

            ## Normalized N-index:

            Nindex <- list(c(1, 4, 1), NULL, 1)
            ## Same as a[c(1, 4, 1), , 1, drop=FALSE]:
            DelayedArray:::subset_by_Nindex(a, Nindex)

            Nindex <- list(integer(0), NULL, 1)
            ## Same as a[integer(0), , 1, drop=FALSE]:
            DelayedArray:::subset_by_Nindex(a, Nindex)

            ## Non-normalized N-index:

            Nindex <- list(-3, NULL, 1)
            Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
            ## Same as a[-3, , 1, drop=FALSE]:
            DelayedArray:::subset_by_Nindex(a, Nindex)

            Nindex <- list(IRanges(2, 4), NULL, 1)
            Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
            ## Same as a[2:4, , 1, drop=FALSE]:
            DelayedArray:::subset_by_Nindex(a, Nindex)

            dimnames(a)[[1]] <- LETTERS[1:4]
            Nindex <- list(c("D", "B"), NULL, 1)
            Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
            ## Same as a[c("D", "B"), , 1, drop=FALSE]:
            DelayedArray:::subset_by_Nindex(a, Nindex)
          }
  }
}

\value{
  \code{Lindex2Mindex} returns an M-index.

  \code{Mindex2Lindex} returns an L-index.
}

\seealso{
  \code{\link[base]{arrayInd}} in the \pkg{base} package.
}

\examples{
dim <- 4:2
Mindex2Lindex(c(4, 3, 1), dim)
Mindex2Lindex(c(4, 3, 2), dim)

Mindex <- rbind(c(1, 1, 1),
                c(2, 1, 1),
                c(3, 1, 1),
                c(4, 1, 1),
                c(1, 2, 1),
                c(1, 1, 2),
                c(4, 3, 2))

Mindex2Lindex(Mindex, dim)

## With a matrix of dimensions:

dims <- rbind(c(4L, 3L),
              c(5L, 3L),
              c(6L, 3L))

Mindex <- rbind(c(1,  2),
                c(1,  2),
                c(1,  2))

Mindex2Lindex(Mindex, dims)

## Sanity checks:

dim <- c(33:30, 45L, 30L)
stopifnot(Mindex2Lindex(rep(1, 6), dim) == 1)
stopifnot(Mindex2Lindex(dim, dim) == prod(dim))

stopifnot(identical(Mindex2Lindex(arrayInd(1:120, 6:4), 6:4), 1:120))
stopifnot(identical(Mindex2Lindex(arrayInd(840:1, 4:7), 4:7), 840:1))
}
\keyword{utilities}
