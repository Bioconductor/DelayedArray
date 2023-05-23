\name{compat}

\alias{t,DelayedArray-method}
\alias{makeNindexFromArrayViewport}
\alias{class:ArrayGrid}
\alias{ArrayGrid-class}
\alias{ArrayGrid}
\alias{DummyArrayGrid}
\alias{RegularArrayGrid}
\alias{ArbitraryArrayGrid}
\alias{extract_array}
\alias{is_sparse}
\alias{read_block}
\alias{write_block}

\title{Functions and classes that have moved to S4Arrays}

\description{
  Some functions and classes that used to be defined in the
  \pkg{DelayedArray} package have been moved to the new \pkg{S4Arrays}
  package in BioC 3.17.
  The corresponding symbols are still exported by the \pkg{DelayedArray}
  package for backward compatibility with existing code.

  WARNING: This is a temporary situation only. Packages that import
  these symbols from \pkg{DelayedArray} must be modified to import them
  from \pkg{S4Arrays} instead.

  These symbols are actually documented in the \pkg{S4Arrays} package.
  See:
  \itemize{
    \item \code{S4Arrays::\link[S4Arrays]{t.Array}}
    \item \code{S4Arrays::\link[S4Arrays]{makeNindexFromArrayViewport}}
    \item \code{S4Arrays::\link[S4Arrays]{ArrayGrid}}
    \item \code{S4Arrays::\link[S4Arrays]{DummyArrayGrid}}
    \item \code{S4Arrays::\link[S4Arrays]{RegularArrayGrid}}
    \item \code{S4Arrays::\link[S4Arrays]{ArbitraryArrayGrid}}
    \item \code{S4Arrays::\link[S4Arrays]{extract_array}}
    \item \code{S4Arrays::\link[S4Arrays]{is_sparse}}
    \item \code{S4Arrays::\link[S4Arrays]{read_block}}
    \item \code{S4Arrays::\link[S4Arrays]{write_block}}
  }
}

\keyword{internal}
