\name{RleArraySeed-class}
\docType{class}

\alias{class:RleArraySeed}
\alias{RleArraySeed-class}
\alias{RleArraySeed}

\alias{class:SolidRleArraySeed}
\alias{SolidRleArraySeed-class}
\alias{SolidRleArraySeed}

\alias{class:RleRealizationSink}
\alias{RleRealizationSink-class}
\alias{RleRealizationSink}

\alias{type,RleRealizationSink-method}

\alias{class:ChunkedRleArraySeed}
\alias{ChunkedRleArraySeed-class}
\alias{ChunkedRleArraySeed}

\alias{dim,RleArraySeed-method}
\alias{dimnames,RleArraySeed-method}

\alias{coerce,SolidRleArraySeed,Rle-method}
\alias{coerce,RleRealizationSink,RleList-method}
\alias{coerce,RleRealizationSink,Rle-method}

\alias{extract_array,SolidRleArraySeed-method}
\alias{extract_array,ChunkedRleArraySeed-method}

\alias{coerce,RleRealizationSink,ChunkedRleArraySeed-method}

\alias{coerce,ChunkedRleArraySeed,SolidRleArraySeed-method}

\title{RleArraySeed objects}

\description{
  RleArraySeed is a low-level helper class for representing an in-memory
  Run Length Encoded array-like dataset.
  RleArraySeed objects are not intended to be used directly. Most
  end users should create and manipulate \link{RleArray} objects
  instead. See \code{?\link{RleArray}} for more information.
}

\details{
  No operation can be performed directly on an RleArraySeed object.
  It first needs to be wrapped in a \link[DelayedArray]{DelayedArray}
  object. The result of this wrapping is an \link{RleArray} object
  (an \link{RleArray} object is just an RleArraySeed object wrapped
  in a \link[DelayedArray]{DelayedArray} object).
}

\seealso{
  \itemize{
    \item \link{RleArray} objects.

    \item \link[S4Vectors]{Rle} objects in the \pkg{S4Vectors} package.
  }
}

\keyword{classes}
\keyword{methods}
