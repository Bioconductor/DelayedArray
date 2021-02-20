\name{ConstantArray}
\docType{class}

\alias{class:ConstantArraySeed}
\alias{ConstantArraySeed-class}
\alias{ConstantArraySeed}

\alias{extract_array,ConstantArraySeed-method}
\alias{extract_sparse_array,ConstantArraySeed-method}
\alias{is_sparse,ConstantArraySeed-method}

\alias{class:ConstantArray}
\alias{ConstantArray-class}
\alias{ConstantArray}

\alias{class:ConstantMatrix}
\alias{ConstantMatrix-class}
\alias{ConstantMatrix}

\alias{matrixClass,ConstantArray-method}
\alias{DelayedArray,ConstantArraySeed-method}
\alias{coerce,ConstantArray,ConstantMatrix-method}
\alias{coerce,ConstantMatrix,ConstantArray-method}

\title{A DelayedArray subclass that contains a constant value}

\description{
  A \link[DelayedArray]{DelayedArray} subclass to efficiently mimic an array
  containing a constant value, without actually creating said array in memory.
}

\usage{
## Constructor function:
ConstantArray(dim, value=NA)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the ConstantArray
    object to create.
  }
  \item{value}{
    Vector (atomic or list) of length 1, containing the value to fill the
    matrix.
  }
}

\value{
  A ConstantArray (or ConstantMatrix) object. (Note that ConstantMatrix
  extends ConstantArray.)
}

\details{
  This class allows us to efficiently create arrays containing a single value.
  For example, we can create matrices full of \code{NA} values, to serve as
  placeholders for missing assays when combining SummarizedExperiment objects.
}

\author{Aaron Lun}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \link{RleArray} objects for representing in-memory Run Length
          Encoded array-like datasets.
  }
}

\examples{
## This would ordinarily take up 8 TB of memory:
CM <- ConstantArray(c(1e6, 1e6), value=NA_real_)
CM

CM2 <-ConstantArray(c(4, 1e6), value=55)
rbind(CM, CM2)
}
\keyword{classes}
\keyword{methods}
