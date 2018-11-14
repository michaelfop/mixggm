\name{penalty}
\alias{penalty}
\title{Penalty functions for graph structure search}

\description{
Collection of penalty functions employed in the structural-EM algorithm and penalized maximum likelihood estimation for graph structure search.
}

\details{
The choice of the penalty function is via argument \code{penalty} in the functions \code{\link{searchGGM}} and \code{\link{mixGGM}}. Possible options are \code{"bic"} (default), \code{"ebic"}, \code{"erdos"}, and \code{"power"}. Functions \code{"ebic"}, \code{"erdos"}, and \code{"power"} depend also on a hyperparameter \code{beta} which can be set using the corresponding argument in \code{\link{searchGGM}} and \code{\link{mixGGM}}.

Let denote with \code{E} the number of non-zero entries in the adjacency matrix corresponding to the graph structure of a covariance/concentration graph model (i.e. the number of edges); \code{N} and \code{V} denote number of observations and number of variables (or nodes). The above options correspond to the following penalty functions:
\itemize{
\item \code{"bic"} -- A BIC-like penalty term is placed on the structure of a graph. This penalty is given by:
\deqn{0.5 E \log N}{0.5*E*log(N)}
The hyperparameter \code{beta} is not used.

\item \code{"ebic"} -- An EBIC-like penalty term for graphical models is placed on the structure of a graph. This penalty is given by:
\deqn{0.5 E \log N + 2\beta E \log V}{0.5*log(N)*E + 2*beta*log(V)*E}
For this penalty function, \code{beta} is a value in the range \code{[0,1]}. Default is \code{beta = 1}, encouraging sparser models. Clearly the case \code{beta = 0} corresponds to \code{"bic"}.

\item \code{"erdos"} -- Let denote by \code{T} the number of all possible edges in a graph, i.e. \eqn{T = {V\choose 2}}{T = choose(V,2)}.  The penalty function is given by:
\deqn{-E\log \beta - (T-E)\log(1 - \beta)}{-E*log(beta) -(T - E)*log(1 - beta)}
For this penalty function, \code{beta} is a value in the range \code{(0,1)}. For small values of \code{beta} the penalization tends to favor situations where the graph decomposes into disjoint blocks. Default is \code{beta = log(V)/T}, a value for which the expected number of arcs is equal to \code{log(V)} and such that the graph will almost surely have disconnected components.

\item \code{"power"} -- Let denote with \code{d_j} the degree of node \code{j}, i. e. the number of nodes connected to it. This penalty function is defined as:
\deqn{\beta \sum_j^V \log(d_j + 1)}{beta * sum_j log(d_j + 1) }
In this case, \code{beta} is a positive value. Default is \code{beta = log(NV)}, a value which place the penalty term on a similar magnitude of \code{"bic"} and \code{"ebic"}, but denser graphs will tend to be less penalized.
}

An user-defined penalty function can be also provided in input of argument \code{penalty} in the functions \code{\link{searchGGM}} and \code{\link{mixGGM}}. In this case, the penalty must be an object of class \code{"function"} and have as argument \code{graph}, like for example \code{"f <- function(graph, beta)"}; see "Examples".

See also \code{\link{searchGGM}} and \code{\link{mixGGM}} for some examples.
}

\examples{

# fit concentration graph model with power law penalty
data(ability.cov)
N <- ability.cov$n.obs
mod1 <- searchGGM(S = ability.cov$cov, N = ability.cov$n.obs,
                  model = "concentration", penalty = "power", beta = 2*log(N))
mod1
plot(mod1)


\dontrun{

# two disconnected blocks of correlated variables
library(MASS)
V <- 10
N <- 500
mu <- rep(0, V)
sigma <- matrix(0.9, V,V)
diag(sigma) <- 1
x <- cbind( MASS::mvrnorm(N, mu, sigma),
            MASS::mvrnorm(N, mu, sigma) )
#
# fit a covariance graph with erdos penalty
mod2 <- searchGGM(x, model = "covariance",
                  penalty = "erdos")
plot(mod2, "adjacency")


# user defined penalty function
data(iris)
x <- iris[,-5]
N <- nrow(x)
V <- ncol(x)
ref <- matrix(0, V, V)
#
# penalize graphs different from a reference graph structure
myPenalty <- function(graph, beta)
{
  beta * sum( abs(graph - ref) )
}
#
mod3 <- mixGGM(x, K = 3, model = "covariance",
              penalty = myPenalty, beta = 2*V*log(N))
plot(mod3)

}

}


\references{
Fop, M., Murphy, T.B., and Scrucca, L. (2018). Model-based clustering with sparse covariance matrices. \emph{Statistics and Computing}. To appear.
}
