\name{SparseData-class}
\docType{class}

\alias{class:SparseData}
\alias{SparseData-class}
\alias{SparseData}

\alias{length,SparseData-method}
\alias{aind}
\alias{aind,SparseData-method}
\alias{nzdata}
\alias{nzdata,SparseData-method}

\alias{dense2sparse}
\alias{sparse2dense}

\alias{as.array.SparseData}
\alias{as.array,SparseData-method}

\alias{as.matrix.SparseData}
\alias{as.matrix,SparseData-method}

\alias{coerce,ANY,SparseData-method}
\alias{coerce,dgCMatrix,SparseData-method}
\alias{coerce,SparseData,dgCMatrix-method}
\alias{coerce,SparseData,sparseMatrix-method}

\title{SparseData objects}

\description{
  SparseData objects are used internally to support block processing
  of array-like objects.
}

\usage{
## Constructor function:
SparseData(dim, aind=NULL, nzdata=NULL, check=TRUE)

## Two low-level utilities:
dense2sparse(x)
sparse2dense(sparse_data)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the
    SparseData object to create.
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
  \item{sparse_data}{
    A SparseData object.
  }
}

\value{
  \itemize{
    \item For \code{SparseData()}: A SparseData instance.

    \item For \code{dense2sparse()}: A SparseData instance.

    \item For \code{sparse2dense()}: An ordinary array.
  }
}

\seealso{
  \itemize{
    \item The \code{\link{read_sparse_block}} function.

    \item \link{block_processing} for more information about block processing
          of an array-like object.

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
sparse_data1 <- SparseData(5:3, aind1, nzdata1)

dim(sparse_data1)
length(sparse_data1)
aind(sparse_data1)
nzdata(sparse_data1)

sparse2dense(sparse_data1)
as.array(sparse_data1)  # same as sparse2dense(sparse_data1)

\dontrun{
as.matrix(sparse_data1)  # error!
}
## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------
m2 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
sparse_data2 <- dense2sparse(m2)
dim(sparse_data2)
length(sparse_data2)
aind(sparse_data2)
nzdata(sparse_data2)

stopifnot(identical(sparse2dense(sparse_data2), m2))

as.matrix(sparse_data2)  # same as sparse2dense(sparse_data2)

## Go back and forth between SparseData and dgCMatrix objects:
M2 <- as(sparse_data2, "dgCMatrix")
stopifnot(identical(M2, as(m2, "dgCMatrix")))

sparse_data2b <- as(M2, "SparseData")
## 'sparse_data2b' is the same as 'sparse_data2' except that
## 'nzdata(sparse_data2b)' is of type numeric instead of integer:
all.equal(sparse_data2b, sparse_data2)
typeof(nzdata(sparse_data2b))  # numeric
typeof(nzdata(sparse_data2))   # integer
}
\keyword{classes}
\keyword{methods}
