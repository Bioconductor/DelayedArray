\name{ConstantArray}
\alias{ConstantMatrix-class}
\alias{ConstantArray}
\alias{ConstantArray-class}
\alias{ConstantArraySeed-class}
\alias{extract_array,ConstantArraySeed-method}
\alias{extract_sparse_array,ConstantArraySeed-method}
\alias{DelayedArray,ConstantArraySeed-method}
\alias{ConstantArraySeed}

\title{A DelayedArray subclass that contains a constant value}

\description{
A \linkS4class{DelayedArray} backend to efficiently mimic a matrix containing a constant value, without actually creating said matrix in memory.
}

\usage{
ConstantArraySeed(dim, value)

ConstantArray(dim, value)

\S4method{DelayedArray}{ConstantArraySeed}(seed)
}

\arguments{
\item{dim}{Integer vector of length greater than 1, containing the dimensions of the array.}

\item{value}{Atomic scalar containing the value to fill the matrix.}

\item{seed}{A ConstantArraySeed object.}
}
\value{
The \code{ConstantArraySeed} constructor returns a ConstantArraySeed object.

The \code{ConstantArray} and \code{DelayedArray} constructors return a ConstantArray object.
}

\details{
This class allows us to efficiently create arrays containing a single value.
For example, we can create matrices full of \code{NA} values, to serve as placeholders for missing assays when combining SummarizedExperiment objects.
We use this class instead of the \linkS4class{RleArray} as the latter requires some workarounds when the product of the dimensions is greater than the maximum integer value.
}

\examples{
# This would ordinarily take up 8 TB of memory:
out <- ConstantArray(c(1e6, 1e6), value=NA_real_)
out
}

\author{
Aaron Lun
}
