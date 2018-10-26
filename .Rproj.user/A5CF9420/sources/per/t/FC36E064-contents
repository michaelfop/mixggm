\name{control-parameters}
\alias{ctrlICF}
\alias{ctrlSTEP}
\alias{ctrlGA}
\alias{ctrlEM}
\alias{ctrlREG}

\title{Set control parameters for various purposes}

\description{
Set control parameters for graphical model estimation, graph structure search via stepwise or genetic algorithm, mixture model fitting via structural-EM algorithm, and Bayesian regularization.
}

\usage{
ctrlICF(tol = 1e-04, maxiter = 1e03)

ctrlSTEP(occamAdd = Inf, occamRem = Inf, start = NULL)

ctrlGA(popSize = 50, pcrossover = 0.8, pmutation = 0.1,
       maxiter = 100, run = maxiter/2,
       elitism = base::max(1, round(popSize*0.05)))

ctrlEM(tol = 1e-05, maxiter = 1e02, subset = NULL, printMsg = FALSE)

ctrlREG(data, K, scaleType = c("full", "fixed", "one", "diag"),
        scale = NULL, psi = NULL)
}


\arguments{
\item{tol}{Tolerance value for judging when convergence has been reached. Used in estimation of a Gaussian graphical model and in the structural-EM algorithm.}

\item{maxiter}{Maximum number of iterations in the Gaussian graphical model estimation algorithm, the genetic algorithm for structure search, and the structural-EM algorithm.}

\item{occamAdd, occamRem}{Set the bounds of the Occam's window for stepwise search. See "Details". Default is \code{Inf}, corresponding to the case of no deletion of candidate graph structures through the search.}

\item{start}{Provide in input a user-defined starting adjacency matrix for stepwise structure search.}

\item{popSize}{Population size. This number corresponds to the number of different graph structures to be considered at each iteration of the genetic algorithm.}

\item{run}{Number of consecutive generations without any improvement in the best fitness value of the structure search procedure before the genetic algorithm is stopped.}

\item{pcrossover}{Probability of crossover between pairs of binary adjacency matrices.}

\item{pmutation}{Probability of mutation in a parent adjacency matrix.}

\item{elitism}{Number of best fitness graph structures to survive at each iteration of the genetic algorithm in the graph structure search procedure.}

\item{subset}{A logical or numeric vector specifying a subset of the data to be used in the initial hierarchical clustering phase employed in the initialization of the structural-EM algorithm. By default no subset is used.}

\item{printMsg}{A logical value indicating whether or not certain warnings (usually related to singularity) should be issued.}

\item{data}{A matrix or data frame of observations. Categorical variables are not allowed. Rows correspond to observations and columns correspond to variables.}

\item{K}{The number of mixture components.}

\item{scaleType}{The type of scale hyperparameter for the prior on the covariance matrix in the case of Bayesian regularization for Gaussian covariance graph model. See "Details". Default is \code{"full"}.}
\item{scale}{The scale hyperparameter for the prior on the covariance matrix in the case of Bayesian regularization for Gaussian covariance graph model.}
\item{psi}{The degrees of freedom hyperparameter for the prior on the covariance matrix in the case of Bayesian regularization for Gaussian covariance graph model.}
}


\details{
Function \code{ctrlICF} is used to set control parameters of the algorithms employed to estimate a Gaussian covariance or concentration graph model.

Function \code{ctrlSTEP} mainly controls the Occam's window used in the stepwise graph structure search. Default is \code{Inf}, corresponding to no Occam's window reduction is implemented. The rationale of the Occam's window is too reduce the space of candidate adjacency matrices during the search by discarding those with a value of the penalized objective function that is too distant from the current optimal value. Graph candidate structures whose penalized objective function value is within \code{occamAdd} from the current optimal value are considered in the next edge-add step of the search. Likewise, graph candidate structures whose penalized objective function value is within \code{occamRem} from the current optimal value are considered in the next edge-remove step of the search. Small values for \code{occamRem} and \code{occamAdd} significantly reduce the space of candidate solutions and the computational cost of the greedy search.

Function \code{ctrlGA} sets parameters of the genetic algorithm used for graph structure search. Arguments correspond to those of function \code{\link[GA]{ga}} in the \code{\link[GA]{GA}} package.

Function \code{ctrlEM} controls standard parameters of the structural-EM algorithm.

Function \code{ctrlREG} is used to set hyperparameters of the conjugate prior on the covariance matrix for Bayesian regularization in Gaussian covariance graph models. The function creates a list of hyperparameters to be given in input in argument \code{regHyperPar} of functions \code{fitGGM}, \code{searchGGM} and \code{mixGGM}.
If not provided in input, the scale hyperparameter is computed on the data via the sample covariance matrix according to the argument \code{scaleType}. If \code{scaleType = "full"}, the scale matrix is proportional to the data covariance matrix; if \code{scaleType = "fixed"}, the scale matrix has determinant equal to \eqn{(\frac{0.001}{K})^{(1/V)}}{(0.001/K)^(1/V)} (Baudry, Celeux, 2015; Fop et al. 2018); if \code{scaleType = "one"} the scale matrix has determinant 1; if \code{scaleType = "diag"} the scale matrix is diagonal. Note that in the case \code{V > N}, if not provided in input the scale matrix is forced to be diagonal.
The hyperparameter \code{psi} controlling the degrees of freedom is set to \code{V + 2} by default.
}

\value{
A list of parameters values.
}


\examples{

\dontrun{

# ga search with increased mutation probability
data(banknote, package = "mclust")
mod1 <- searchGGM(banknote[,-1], model = "concentration", search = "ga",
                  ctrlGa = ctrlGA(pmutation = 0.3))


# regularization
library(MASS)
V <- 10
N <- 20
mu <- rep(0, V)
sigma <- matrix(0.9, V,V)
diag(sigma) <- 1
x <- cbind( MASS::mvrnorm(N, mu, sigma),
            MASS::mvrnorm(N, mu, sigma),
            MASS::mvrnorm(N, mu, sigma))  # high-dimensional data V = 30, N = 20
#
hyperPar <- ctrlREG(x, K = 1, scaleType = "diag")
mod2 <- searchGGM(x, model = "covariance", penalty = "ebic")   # throws an error
mod2 <- searchGGM(x, model = "covariance", penalty = "ebic",   # regularization
                  regularize = TRUE, regHyperPar = hyperPar)
plot(mod2, "adjacency")


# occam's window
library(MASS)
V <- 20
N <- 500
mu <- rep(0, V)
sigma <- matrix(0.9, V,V)
diag(sigma) <- 1
edges <- rbinom(choose(V,2), 1, 0.3)
A <- matrix(0, V,V)
A[lower.tri(A)] <- edges
A <- A + t(A)
fit <- fitGGM(S = sigma, N = N, graph = A, model = "concentration",
              ctrlIcf = ctrlICF(tol = 1e-06))
sigma <- fit$sigma
#
x <- MASS:::mvrnorm(N, mu, sigma)
#
mod3 <- searchGGM(x, model = "concentration", search = "step-back",
                  ctrlStep = ctrlSTEP(occamAdd = 5, occamRem = 5))
par(mfrow = c(1,2))
plot(fit, what = "adjacency")
plot(mod3, what = "adjacency")

}

}


\references{
Baudry, J.P. and Celeux, G. (2015). EM for mixtures: Initialization requires special care. \emph{Statistics and Computing}, 25(4):713-726.

Fop, M., Murphy, T.B., and Scrucca, L. (2018). Model-based clustering with sparse covariance matrices. \emph{Statistics and Computing}. To appear.
}
