\name{searchGGM}
\alias{searchGGM}
\title{Structure search and estimation for Gaussian graphical models}

\description{
Graph structure search and estimation for Gaussian covariance and concentration graph models.
}

\usage{
searchGGM(data = NULL,
          S = NULL, N = NULL,
          model = c("covariance", "concentration"),
          search = c("step-forw", "step-back", "ga"),
          penalty = c("bic", "ebic", "erdos", "power"),
          beta = NULL,
          start = NULL,
          regularize = FALSE, regHyperPar = NULL,
          ctrlStep = ctrlSTEP(), ctrlGa = ctrlGA(), ctrlIcf = ctrlICF(),
          parallel = FALSE,
          verbose = FALSE, ...)
}

\arguments{
\item{data}{A dataframe or matrix, where rows correspond to observations and columns to variables. Categorical variables are not allowed.}

\item{S}{The sample covariance matrix of the data. If \code{S = NULL}, the maximum likelihood estimate of the covariance matrix is used in the estimation of the graphical model.}

\item{N}{The number of observations. If \code{data = NULL} and \code{S} is provided in input, \code{N} must be provided in input as well.}

\item{model}{The type of Gaussian graphical model. Default is \code{"covariance"}. See "Details".}

\item{search}{The type of structure search algorithm. If \code{search = "step-forw"}, a greedy forward-stepwise search is used to find the optimal graph association structure. If \code{search = "step-back"}, a greedy backward-stepwise search is implemented. If \code{search = "ga"} a stochastic search based on a genetic algorithm is employed. Default is \code{"step-forw"}.}

\item{penalty}{The penalty function used to define a criterion for scoring the candidate graph configurations. Default is \code{"bic"}. See "Details" and \code{\link{penalty}}.}

\item{beta}{The hyperparameter of the penalty function. See "Details" and \code{\link{penalty}}.}

\item{start}{A starting matrix for the estimation algorithm. If \code{NULL}, the starting value is the diagonal sample covariance matrix. Used only when \code{model = "covariance"}.}

\item{regularize}{A logical argument indicating if Bayesian regularization should be performed. Default to \code{FALSE}. Used only when \code{model = "covariance"}.}

\item{regHyperPar}{A list of hyper parameters for Bayesian regularization. Only used when\cr \code{regularization = TRUE}; see also \code{\link{ctrlREG}}.}

\item{ctrlStep}{A list of control parameters used in the stepwise search; see also \code{\link{ctrlSTEP}}.}
\item{ctrlGa}{A list of control parameters for the genetic algorithm; see also \code{\link{ctrlGA}}.}
\item{ctrlIcf}{A list of control parameters employed in the algorithm for estimation of graphical model parameters; see also \code{\link{ctrlICF}}.}

\item{parallel}{A logical argument indicating if parallel computation should be used for structure search. If TRUE, all the available cores are used. The argument could also be set to a numeric integer value specifying the number of cores to be employed.}

\item{verbose}{A logical argument controlling whether iterations of the structure searching and estimation procedure need to be shown or not.}

\item{...}{Additional internal arguments not to be provided by the user.}
}


\value{
An object of class \code{'fitGGM'} containing the optimal estimated marginal or conditional independence Gaussian graphical model.

The output is a list containing:
\item{sigma}{The estimated covariance matrix.}
\item{omega}{The estimated concentration (inverse covariance) matrix.}
\item{graph}{The adjacency matrix corresponding to the optimal marginal or conditional independence graph.}
\item{model}{Estimated model type, whether \code{"covariance"} or \code{"concentration"}.}
\item{loglikPen}{Value of the maximized penalized log-likelihood.}
\item{loglik}{Value of the maximized log-likelihood.}
\item{nPar}{Number of estimated parameters.}
\item{N}{Number of observations.}
\item{V}{Number of variables, corresponding to the number of nodes in the graph.}
\item{penalty}{The type of penalty on the graph structure.}
\item{search}{The search method used for graph structure search.}
\item{GA}{An object of class \code{'ga-class'} with information about the genetic algorithm. Only present when \code{search = "ga"}. See \code{\link[GA]{ga}}.}
}


\details{
The function performs graph association structure search and maximum penalized likelihood estimation of the optimal Gaussian graphical model given the data provided in input.

A Gaussian covariance graph model is estimated if \code{model = "covariance"}, while estimation of a Gaussian covariance graph model is performed if \code{model = "concentration"}. A Gaussian covariance graph model postulates that some variables are marginally independent according to the inferred graph structure. On the other hand, in a Gaussian concentration graph model, variables are conditionally independent given their neighbors in the inferred graph. See also \code{\link{fitGGM}}.

Search for the optimal graph structure and parameter estimation is carried out by maximization of a Gaussian penalized likelihood, given as follows:
\deqn{ \mbox{Covariance:}~~~ \mbox{argmax}_{\Sigma,A} ~~ \ell(X | \Sigma, A) - P(A, \beta) \quad \Sigma \in C^{+}_G(A)}{Covariance:  argmax_(Sigma, A)  L(X | Sigma, A) - P_beta(A) with Sigma in C_G(A) }
\deqn{ \mbox{Concentration:}~~~ \mbox{argmax}_{\Omega,A} ~~ \ell(X | \Omega, A) - P(A, \beta) \quad \Omega \in C^{+}_G(A)}{Concentration:  argmax_(Omega, A)  L(X | Omega, A) - P_beta(A) with Omega in C_G(A) }
where \eqn{C^{+}_G(A)}{C_G(A)} is the collection of sparse positive definite matrices whose zero patterns are given by graph \eqn{G}{G} represented by the adjacency matrix \eqn{A}{A}.

The penalty function \eqn{P(A, \beta)}{P_beta(A)} depends on the structure of graph \eqn{G}{G} through the adjacency matrix \eqn{A}{A} and a parameter \eqn{\beta}{beta}; see \code{\link{penalty}} on how to specify the penalization term and for further information.

For this type of penalized log-likelihood, graph structure search and parameter estimation is a maximization combinatorial problem. For a given candidate structure (i.e. adjacency matrix), association parameters in the covariance or concentration matrix are estimated using the estimation algorithms implemented in \code{\link{fitGGM}}. Regarding structure search, this can be carried out either using a greedy forward-stepwise or a greedy backward-stepwise algorithm, by setting \code{search = "step-forw"} or \code{search = "step-back"} respectively. Alternatively, a stochastic search via genetic algorithm can be used by setting \code{search = "ga"}. The procedure for the forward stepwise search is described in Fop et al. (2018), and the backward is implemented in a similar way; the genetic algorithm procedure relies on the \code{\link[GA]{GA}} package. All the structure searching methods can be run in parallel on a multi-core machine by setting the argument \code{parallel = TRUE}.
}

\examples{

# fit covariance graph model with default forward-stepwise search
data(mtcars)
x <- mtcars[,c(1,3:7)]
mod1 <- searchGGM(x, model = "covariance")
mod1
plot(mod1)
#
# prefer a sparser model
mod2 <- searchGGM(x, model = "covariance", penalty = "ebic")
mod2
plot(mod2)


# fit concentration graph model with backward-stepwise structure search
# with a covariance matrix in input
data(ability.cov)
mod3 <- searchGGM(S = ability.cov$cov, N = ability.cov$n.obs,
                  model = "concentration", search = "step-back")
mod3
mod3$graph
mod3$omega
plot(mod3)

\dontrun{

# generate data from a Markov model
N <- 1000
V <- 20
dat <- matrix(NA, N, V)
dat[,1] <- rnorm(N)
for ( j in 2:V ) dat[,j] <- dat[,j-1] + rnorm(N, sd = 0.5)
mod4 <- searchGGM(data = dat, model = "concentration")        # recover the model
plot(mod4, what = "adjacency")

}
}


\references{
Fop, M., Murphy, T.B., and Scrucca, L. (2018). Model-based clustering with sparse covariance matrices. \emph{Statistics and Computing}. To appear.

Scrucca, L. (2017). On some extensions to GA package: Hybrid optimisation, parallelisation and islands evolution. \emph{The R Journal}, 9(1), 187-206.

Scrucca, L. (2013). GA: A package for genetic algorithms in R. \emph{Journal of Statistical Software}, 53(4), 1-3.
}
