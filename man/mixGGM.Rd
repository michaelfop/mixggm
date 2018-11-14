\name{mixGGM}
\alias{mixGGM}
\alias{print.mixGGM}
\title{Mixture of Gaussian graphical models}

\description{
Estimation of a mixture of Gaussian covariance or concentration graph models using structural-EM algorithm. The mixture model returned is the optimal model according to BIC.
}

\usage{
mixGGM(data, K = 1:3,
       model = c("covariance", "concentration"),
       search = c("step-forw", "step-back", "ga"),
       penalty = c("bic", "ebic", "erdos", "power"),
       beta = NULL,
       regularize = FALSE, regHyperPar = NULL,
       ctrlEm = ctrlEM(),
       ctrlStep = ctrlSTEP(), ctrlGa = ctrlGA(),
       ctrlIcf = ctrlICF(),
       keepAll = FALSE,
       parallel = FALSE,
       verbose = TRUE)}

\arguments{
\item{data}{A dataframe or matrix, where rows correspond to observations and columns to variables. Categorical variables are not allowed.}

\item{K}{An integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated.}

\item{model}{The type of Gaussian graphical model. Default is \code{"covariance"}. See "Details".}

\item{search}{The type of structure search algorithm. If \code{search = "step-forw"}, a greedy forward-stepwise search is used to find the optimal graph association structure. If \code{search = "step-back"}, a greedy backward-stepwise search is implemented. If \code{search = "ga"} a stochastic search based on a genetic algorithm is employed. Default is \code{"step-forw"}.}

\item{penalty}{The penalty function used to define a criterion for scoring the candidate graph configurations. Default is \code{"bic"}. See "Details" and \code{\link{penalty}}.}

\item{beta}{The hyperparameter of the penalty function. See "Details" and \code{\link{penalty}}.}

\item{regularize}{A logical argument indicating if Bayesian regularization should be performed. Default to \code{FALSE}. Used only when \code{model = "covariance"}.}

\item{regHyperPar}{A list of hyper parameters for Bayesian regularization. Only used when\cr \code{regularization = TRUE}; see also \code{\link{ctrlREG}}.}

\item{ctrlEm}{A list of control parameters used in the structural-EM algorithm; see also \code{\link{ctrlEM}}.}
\item{ctrlStep}{A list of control parameters used in the stepwise search; see also \code{\link{ctrlSTEP}}.}
\item{ctrlGa}{A list of control parameters for the genetic algorithm; see also \code{\link{ctrlGA}}.}
\item{ctrlIcf}{A list of control parameters employed in the algorithm for estimation of graphical model parameters; see also \code{\link{ctrlICF}}.}

\item{keepAll}{A logical argument. If \code{TRUE}, also all the mixture models estimated for the values of \code{K} given in input are returned}

\item{parallel}{A logical argument indicating if parallel computation should be used for structure search in the M step of the structural-EM algorithm. If TRUE, all the available cores are used. The argument could also be set to a numeric integer value specifying the number of cores to be employed.}

\item{verbose}{If \code{TRUE} a progress bar will be shown.}
}


\value{
An object of class \code{'mixGGM'} containing the optimal estimated mixture of Gaussian graphical models.

The output is a list containing:
\item{parameters}{A list with the following components:
\describe{
        \item{\code{tau}}{A vector containing the estimated mixing proportions.}
        \item{\code{mu}}{The mean for each mixture component. Columns denote the mixture components.}
        \item{\code{sigma}}{An array containing the components covariance matrices.}
        \item{\code{omega}}{An array containing the components concentration (inverse covariance) matrices.}
  }
}
\item{graph}{An array with the adjacency matrices corresponding to the optimal marginal or conditional independence graphs for each mixture component.}
\item{N}{Number of observations in the data.}
\item{V}{Number of variables in the data, corresponding to the number of nodes in the graphs.}
\item{K}{Number of selected mixture components.}
\item{loglik}{Value of the maximized log-likelihood.}
\item{loglikPen}{Value of the maximized penalized log-likelihood.}
\item{loglikReg}{Value of the maximized regularized log-likelihood. If \code{regularize = FALSE}, this value is equal to \code{loglik}}
\item{nPar}{A vector with two entries:
\describe{
        \item{\code{depPar}}{Total number of dependence parameters. If \code{model =} \code{"covariance"}, this is the total number of non-zero covariance parameters, while if \code{model =} \code{"concentration"}, it corresponds to the total number of non-zero concentration parameters.}
        \item{\code{totPar}}{Total number of mixture parameters.}
  }
}
\item{z}{A matrix whose \code{[i,k]}th entry is the probability that observation \code{i} of the data belongs to the \code{k}th class.}
\item{classification}{Classification corresponding to the maximum a posteriori of matrix \code{z}.}
\item{bic}{Optimal BIC value.}
\item{BIC}{All BIC values.}
\item{data}{The data matrix provided in input.}
\item{model}{Estimated model type, whether \code{"covariance"} or \code{"concentration"}.}
\item{penalty}{The type of penalty on the graph structure.}
\item{search}{The search method used for graph structure search.}
\item{keepAll}{A list containing all the estimated models. Provided in output only when \code{keepAll = TRUE}.}
}


\details{
Estimation of a mixture of Gaussian graphical models by means of maximization of a penalized log-likelihood via structural-EM algorithm. The mixture model in output is the optimal model selected by BIC.

If \code{model = "covariance"}, a \strong{mixture of Gaussian covariance graph models} is estimated. The Gaussian mixture is parameterized in terms of the components covariance matrices and the component adjacency matrices correspond to \emph{marginal independence} constraints among the variables:
\deqn{X \sim \sum_{k}^K\tau_k\mathcal{N}(\mu_k, \Sigma_k) \quad \Sigma_k \in C^{+}_G(A_k)}{X ~ sum_k tau_k N(mu_k, Sigma_k) with Sigma_k in C_G(A_k) }
Variables associated to two non-connected edges in the graphs are marginally independent and have different marginal association patterns across the mixture components. As a result, the covariance matrices \code{sigma} are estimated to be sparse according to the inferred graph structures.

If \code{model = "concentration"}, estimation of a \strong{mixture of Gaussian concentration graph model} is performed. The Gaussian mixture is parameterized in terms of the components concentration matrices and the component adjacency matrices correspond to \emph{conditional independence} constraints among the variables:
\deqn{X \sim \sum_{k}^K\tau_k\mathcal{N}(\mu_k, \Omega_k) \quad \Omega_k \in C^{+}_G(A_k)}{X ~ sum_k N(mu_k, Omega_k) with Omega_k in C_G(A_k) }
Variables associated to two non-adjacent edges in the graph are conditionally independent given their common neighbors and have different conditional dependence patterns across the mixture components. It results in the concentration matrices \code{omega} being estimated to be sparse according to the inferred graph structures.

Arguments \code{penalty} and \code{search} are used to define the type of penalty on the graph configuration and the structure search method in the structural-EM algorithm. The penalization term depends on the hyperparameter \code{beta} according to the type of penalty function. See \code{\link{searchGGM}} and \code{\link{penalty}} for more details.
}

\examples{

# fit a mixture of concentration graph models
data(iris)
mod1 <- mixGGM(iris[,-5], model = "concentration")
plot(mod1, what = "graph")
plot(mod1, what = "classification")


\dontrun{

# a simple simulated data example
library(MASS)
N <- 200
tau <- c(0.3, 0.7)
Nk <- rowSums( rmultinom(N, 1, tau) )
class <- rep(1:2, Nk)
sigma1 <- diag(2)                          # independent variables
sigma2 <- matrix( c(1,0.9,0.9,1), 2,2 )    # correlated variables
mu1 <- c(0, 0)
mu2 <- c(2, 3)
x <- rbind( MASS::mvrnorm(Nk[1], mu1, sigma1),
            MASS::mvrnorm(Nk[2], mu2, sigma2)
)
mod2 <- mixGGM(x)
plot(mod2)
plot(mod2, what = "classification")

# fit a mixture of covariance graph models
data(wine, package = "gclus")
mod3 <- mixGGM(wine[,-1], K = 1:4, model = "covariance",
               penalty = "erdos", beta = 0.01)
plot(mod3, what = "graph")
plot(mod3, what = "classification", dimens = 1:4)


# complex simulated data example
N <- 500
V <- 20
tau <- c(0.3, 0.7)
Nk <- rowSums( rmultinom(N, 1, tau) )
class <- rep(1:2, Nk)
sigma1 <- rWishart(1, V+1, diag(V))[,,1]
mu1 <- rep(0, V)
mu2 <- rnorm(V, 0.5, 2)
x1 <- MASS::mvrnorm(Nk[1], mu1, sigma1)
x2 <- matrix(NA, Nk[2], V)
x2[,1] <- rnorm(Nk[2])
for ( j in 2:V ) x2[,j] <- x2[,j-1] + rnorm(Nk[2], mu2[j], sd = 0.5)
x <- rbind(x1, x2)
#
mod4 <- mixGGM(x, K = 1:4, model = "concentration",
               penalty = "ebic", beta = 0.5)
plot(mod4, what = "classification", dimens = c(1,5,10,15,20) )
plot(mod4, what = "graph")
plot(mod4, what = "adjacency")
table(class, mod4$classification)
#
mc <- mclust::Mclust(x, G = 1:4)
mc$bic
mod4$bic

}

}


\references{
Fop, M., Murphy, T.B., and Scrucca, L. (2018). Model-based clustering with sparse covariance matrices. \emph{Statistics and Computing}. To appear.
}
