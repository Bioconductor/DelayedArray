\name{SparseArraySeed-class}
\docType{class}

\alias{class:SparseArraySeed}
\alias{SparseArraySeed-class}
\alias{SparseArraySeed}

\alias{aind}
\alias{aind,SparseArraySeed-method}
\alias{nzdata}
\alias{nzdata,SparseArraySeed-method}
\alias{sparsity}
\alias{sparsity,SparseArraySeed-method}

\alias{dense2sparse}
\alias{sparse2dense}

\alias{is_sparse,ANY-method}
\alias{extract_sparse_array}

\alias{is_sparse,SparseArraySeed-method}
\alias{extract_sparse_array,SparseArraySeed-method}
\alias{extract_array,SparseArraySeed-method}

\alias{is_sparse,dgCMatrix-method}
\alias{extract_sparse_array,dgCMatrix-method}

\alias{as.array.SparseArraySeed}
\alias{as.array,SparseArraySeed-method}

\alias{as.matrix.SparseArraySeed}
\alias{as.matrix,SparseArraySeed-method}

\alias{coerce,ANY,SparseArraySeed-method}

\alias{coerce,dgCMatrix,SparseArraySeed-method}
\alias{coerce,SparseArraySeed,dgCMatrix-method}
\alias{coerce,SparseArraySeed,sparseMatrix-method}

\alias{aperm.SparseArraySeed}
\alias{aperm,SparseArraySeed-method}

\title{SparseArraySeed objects}

\description{
  SparseArraySeed objects are used internally to support block processing
  of array-like objects.
}

\usage{
## Constructor function:
SparseArraySeed(dim, aind=NULL, nzdata=NULL, check=TRUE)

## Getters (in addition to dim() and length()):
aind(x)
nzdata(x)
sparsity(x)

## Two low-level utilities:
dense2sparse(x)
sparse2dense(sas)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the
    SparseArraySeed object to create.
  }
  \item{aind}{
    A matrix containing the array indices of the nonzero data.

    This must be an integer matrix like one returned by
    \code{base::\link[base]{arrayInd}}, that is, with \code{length(dim)}
    columns and where each row is an n-uplet representing an \emph{array index}.
  }
  \item{nzdata}{
    A vector of length \code{nrow(aind)} containing the nonzero data.
  }
  \item{check}{
    Should the object be validated upon construction?
  }
  \item{x}{
    A SparseArraySeed object for the \code{aind}, \code{nzdata}, and
    \code{sparsity} getters.

    An array-like object for \code{dense2sparse}.
  }
  \item{sas}{
    A SparseArraySeed object.
  }
}

\value{
  \itemize{
    \item For \code{SparseArraySeed()}: A SparseArraySeed instance.

    \item For \code{aind()}: The matrix containing the array indices of the
          nonzero data.

    \item For \code{nzdata()}: The vector of nonzero data.

    \item For \code{sparsity()}: The number of zero-valued elements
          in the implicit array divided by the total number of array
          elements (a.k.a. the length of the array).

    \item For \code{dense2sparse()}: A SparseArraySeed instance.

    \item For \code{sparse2dense()}: An ordinary array.
  }
}

\seealso{
  \itemize{
    \item The \code{\link{read_sparse_block}} function.

    \item \code{\link{blockApply}} for more information about block processing
          of array-like objects.

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
sas1 <- SparseArraySeed(5:3, aind1, nzdata1)

dim(sas1)        # the dimensions of the implicit array
length(sas1)     # the length of the implicit array
aind(sas1)
nzdata(sas1)
sparsity(sas1)

sparse2dense(sas1)
as.array(sas1)   # same as sparse2dense(sas1)

\dontrun{
as.matrix(sas1)  # error!
}
## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------
m2 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
sas2 <- dense2sparse(m2)
dim(sas2)
length(sas2)
aind(sas2)
nzdata(sas2)
sparsity(sas2)

stopifnot(identical(sparse2dense(sas2), m2))

as.matrix(sas2)  # same as sparse2dense(sas2)

t(sas2)
stopifnot(identical(as.matrix(t(sas2)), t(as.matrix(sas2))))

## Go back and forth between SparseArraySeed and dgCMatrix objects:
M2 <- as(sas2, "dgCMatrix")
stopifnot(identical(M2, as(m2, "dgCMatrix")))

sas2b <- as(M2, "SparseArraySeed")
## 'sas2b' is the same as 'sas2' except that
## 'nzdata(sas2b)' is of type numeric instead of integer:
all.equal(sas2b, sas2)
typeof(nzdata(sas2b))  # numeric
typeof(nzdata(sas2))   # integer

## ---------------------------------------------------------------------
## SEED CONTRACT
## ---------------------------------------------------------------------
## SparseArraySeed objects comply with the "seed contract".
## In particular they support extract_array():
extract_array(sas1, list(c(5, 3:2, 5), NULL, 3))

## See '?extract_array' for more information about the "seed contract".

## This means that they can be wrapped in a DelayedArray object:
A1 <- DelayedArray(sas1)
A1

## A big very sparse DelayedMatrix object:
aind3 <- cbind(sample(25000, 600000, replace=TRUE),
               sample(195000, 600000, replace=TRUE))
nzdata3 <- runif(600000)
sas3 <- SparseArraySeed(c(25000, 195000), aind3, nzdata3)
sparsity(sas3)

M3 <- DelayedArray(sas3)
M3
colSums(M3[ , 1:20])
}
\keyword{classes}
\keyword{methods}
