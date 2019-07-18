\name{summary.mixGGM}
\alias{print.mixGGM}
\alias{summary.mixGGM}
\alias{print.summary.mixGGM}

\title{Summarizing Mixture of Gaussian Graphical Models}

\description{Summary method for class \code{'mixGGM'}.}

\usage{
\method{summary}{mixGGM}(object, graphs = TRUE, clusters = FALSE, parameters = FALSE, \dots)
\method{print}{summary.mixGGM}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{'mixGGM'} resulting of a call to \code{\link{mixGGM}}.}

  \item{x}{An object of class \code{'summary.mixGGM'}, usually, a result of a call to \code{summary.mixGGM}.}

  \item{graphs}{Logical; if \code{TRUE}, the adjancecy matrices of mixture graphs is printed.}

  \item{clusters}{Logical; if \code{TRUE}, the table of MAP classification/clustering of observations is printed.}

  \item{parameters}{Logical; if \code{TRUE}, the parameters of mixture components are printed.}

  \item{digits}{The number of significant digits to use when printing.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

% \details{}

% \value{}

% \note{}

\seealso{\code{\link{mixGGM}}.}

% \examples{}
