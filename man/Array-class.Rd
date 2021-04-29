\name{Array-class}
\docType{class}

\alias{class:Array}
\alias{Array-class}
\alias{Array}

\alias{length,Array-method}
\alias{isEmpty,Array-method}
\alias{[[,Array-method}

\alias{t.Array}
\alias{t,Array-method}

\alias{coerce,Array,dgCMatrix-method}
\alias{coerce,Array,dgRMatrix-method}
\alias{coerce,Array,lgCMatrix-method}
\alias{coerce,Array,lgRMatrix-method}
\alias{coerce,Array,CsparseMatrix-method}
\alias{coerce,Array,RsparseMatrix-method}
\alias{coerce,Array,sparseMatrix-method}

\alias{showAsCell,Array-method}

\title{Array objects}

\description{
  Array is a virtual class intended to be extended by concrete subclasses
  with an array-like semantic.
}

\seealso{
  \link{DelayedArray}, \link{ArrayGrid}, and \link{ArrayViewport} for
  examples of classes with an array-like semantic.
}

\examples{
showClass("Array")  # virtual class with no slots
}
\keyword{internal}
