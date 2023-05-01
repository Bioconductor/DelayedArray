\name{compat}

\alias{t,DelayedArray-method}
\alias{makeNindexFromArrayViewport}
\alias{DummyArrayGrid}
\alias{RegularArrayGrid}
\alias{ArbitraryArrayGrid}
\alias{extract_array}
\alias{is_sparse}
\alias{read_block}
\alias{write_block}

\title{Symbols defined in S4Arrays and re-exported by DelayedArray}

\description{
  Some symbols defined in the \pkg{S4Arrays} package are currently
  re-exported by the \pkg{DelayedArray} package.

  WARNING: This is a temporary situation only. Packages that import
  these symbols from \pkg{DelayedArray} must be modified to import them
  from \pkg{S4Arrays} instead.

  These symbols are actually documented in the \pkg{S4Arrays} package.
  See:
  \itemize{
    \item \code{S4Arrays::\link[S4Arrays]{t.Array}}
    \item \code{S4Arrays::\link[S4Arrays]{makeNindexFromArrayViewport}}
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
