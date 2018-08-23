\name{SparseArray-class}
\docType{class}

\alias{class:SparseArray}
\alias{SparseArray-class}
\alias{SparseArray}

\alias{aind}
\alias{aind,SparseArray-method}
\alias{nzdata}
\alias{nzdata,SparseArray-method}

\alias{dense2sparse}
\alias{sparse2dense}

\alias{extract_array,SparseArray-method}

\alias{as.array.SparseArray}
\alias{as.array,SparseArray-method}

\alias{as.matrix.SparseArray}
\alias{as.matrix,SparseArray-method}

\alias{coerce,ANY,SparseArray-method}

\alias{coerce,dgCMatrix,SparseArray-method}
\alias{coerce,SparseArray,dgCMatrix-method}
\alias{coerce,SparseArray,sparseMatrix-method}

\title{SparseArray objects}

\description{
  SparseArray objects are used internally to support block processing
  of array-like objects.
}

\usage{
## Constructor function:
SparseArray(dim, aind=NULL, nzdata=NULL, check=TRUE)

## Two low-level utilities:
dense2sparse(x)
sparse2dense(sparse_array)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the
    SparseArray object to create.
  }
  \item{aind}{
    An integer matrix like one returned by \code{base::\link[base]{arrayInd}},
    that is, with \code{length(dim)} columns and where each row is an
    n-uplet representing an array index.
  }
  \item{nzdata}{
    A vector of length \code{nrow(aind)} containing the nonzero data.
  }
  \item{check}{
    Should the object be validated upon construction?
  }
  \item{x}{
    An array-like object.
  }
  \item{sparse_array}{
    A SparseArray object.
  }
}

\value{
  \itemize{
    \item For \code{SparseArray()}: A SparseArray instance.

    \item For \code{dense2sparse()}: A SparseArray instance.

    \item For \code{sparse2dense()}: An ordinary array.
  }
}

\seealso{
  \itemize{
    \item The \code{\link{read_sparse_block}} function.

    \item \link{block_processing} for more information about block processing
          of an array-like object.

    \item \code{\link{extract_array}}.

    \item \link{DelayedArray} objects.

    \item \code{\link[base]{arrayInd}} in the \pkg{base} package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## EXAMPLE 1
## ---------------------------------------------------------------------
aind1 <- rbind(c(2,4,3), c(2,1,3), c(5,4,3), c(5,3,3),
               c(5,4,1), c(5,1,1), c(5,4,2), c(5,4,1))
nzdata1 <- 11.11 * seq_len(nrow(aind1))
sparse_array1 <- SparseArray(5:3, aind1, nzdata1)

dim(sparse_array1)
length(sparse_array1)
aind(sparse_array1)
nzdata(sparse_array1)

sparse2dense(sparse_array1)
as.array(sparse_array1)   # same as sparse2dense(sparse_array1)

\dontrun{
as.matrix(sparse_array1)  # error!
}
## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------
m2 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
sparse_array2 <- dense2sparse(m2)
dim(sparse_array2)
length(sparse_array2)
aind(sparse_array2)
nzdata(sparse_array2)

stopifnot(identical(sparse2dense(sparse_array2), m2))

as.matrix(sparse_array2)  # same as sparse2dense(sparse_array2)

## Go back and forth between SparseArray and dgCMatrix objects:
M2 <- as(sparse_array2, "dgCMatrix")
stopifnot(identical(M2, as(m2, "dgCMatrix")))

sparse_array2b <- as(M2, "SparseArray")
## 'sparse_array2b' is the same as 'sparse_array2' except that
## 'nzdata(sparse_array2b)' is of type numeric instead of integer:
all.equal(sparse_array2b, sparse_array2)
typeof(nzdata(sparse_array2b))  # numeric
typeof(nzdata(sparse_array2))   # integer

## ---------------------------------------------------------------------
## SEED CONTRACT
## ---------------------------------------------------------------------
## SparseArray comply with the "seed contract". In particular they
## support extract_array():
extract_array(sparse_array1, list(c(5, 3:2, 5), NULL, 3))

## See '?extract_array' for more information about the "seed contract".

## This means that they can be wrapped in a DelayedArray object:
A1 <- DelayedArray(sparse_array1)
A1

## A big very sparse DelayedMatrix object:
aind3 <- cbind(sample(25000, 600000, replace=TRUE),
               sample(195000, 600000, replace=TRUE))
nzdata3 <- runif(600000)
sparse_array3 <- SparseArray(c(25000, 195000), aind3, nzdata3)
M3 <- DelayedArray(sparse_array3)
M3
colSums(M3[ , 1:20])
}
\keyword{classes}
\keyword{methods}
