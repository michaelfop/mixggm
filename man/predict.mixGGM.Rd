\name{predict.mixGGM}
\alias{predict.mixGGM}

\title{Cluster prediction by Mixture of Gaussian Graphical Models}

\description{Cluster prediction for multivariate observations based on Mixture of Gaussian Graphical Models estimated by \code{\link{mixGGM}}.}

\usage{
\method{predict}{mixGGM}(object, newdata, \dots)
}

\arguments{
  \item{object}{An object of class \code{'mixGGM'} resulting from a call to \code{\link{mixGGM}}.}

  \item{newdata}{A data frame or matrix giving the data. If missing the clustering data obtained from the call to \code{\link{mixGGM}} are classified.}

  \item{\dots}{Further arguments passed to or from other methods.}

}

% \details{}

\value{
Returns a list of with the following components:
  \item{classification}{a factor of predicted cluster labels for \code{newdata}.}
  \item{z}{a matrix whose \emph{[i,k]}th entry is the probability that 
           observation \emph{i} in \code{newdata} belongs to the \emph{k}th cluster.}
}

% \note{}

\seealso{\code{\link{mixGGM}}.}

% \examples{}