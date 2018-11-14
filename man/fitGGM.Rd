\name{fitGGM}
\alias{fitGGM}
\alias{print.fitGGM}
\title{Fit a Gaussian graphical model}

\description{
Estimation of a Gaussian graphical model given the graph structure corresponding to marginal or conditional independence restrictions.
}

\usage{
fitGGM(data = NULL,
      S = NULL, N = NULL,
      graph,
      model = c("covariance", "concentration"),
      start = NULL,
      ctrlIcf = ctrlICF(),
      regularize = FALSE,
      regHyperPar = NULL,
      verbose = FALSE, ...)
}

\arguments{
\item{data}{A dataframe or matrix, where rows correspond to observations and columns to variables. Categorical variables are not allowed.}

\item{S}{The sample covariance matrix of the data. If \code{S = NULL}, the maximum likelihood estimate of the covariance matrix is used in the estimation of the graphical model.}

\item{N}{The number of observations. If \code{data = NULL} and \code{S} is provided in input, \code{N} must be provided in input as well.}

\item{graph}{A square symmetric binary adjacency matrix corresponding to the association structure of the graph. See "Details".}

\item{model}{The type of Gaussian graphical model. Default is \code{"covariance"}. See "Details".}

\item{start}{A starting matrix for the estimation algorithm. If \code{NULL}, the starting value is the diagonal sample covariance matrix. Used only when \code{model = "covariance"}.}

\item{ctrlIcf}{A list of control parameters for the numerical algorithm for estimation of graphical model parameters; see also \code{\link{ctrlICF}}.}

\item{regularize}{A logical argument indicating if Bayesian regularization should be performed. Default to \code{FALSE}. Used only when \code{model = "covariance"}.}

\item{regHyperPar}{A list of hyper parameters for Bayesian regularization. Only used when\cr \code{regularization = TRUE}; see also \code{\link{ctrlREG}}.}

\item{verbose}{A logical argument controlling whether iterations of the estimation procedure need to be shown or not.}

\item{...}{Additional internal arguments not to be provided by the user.}
}

\value{
An object of class \code{'fitGGM'} containing the estimated Gaussian graphical model.

The output is a list containing:
\item{sigma}{The estimated covariance matrix.}
\item{omega}{The estimated concentration (inverse covariance) matrix.}
\item{graph}{The adjacency matrix given in input corresponding to the marginal or conditional independence graph.}
\item{model}{Estimated model type, whether \code{"covariance"} or \code{"concentration"}.}
\item{loglik}{Value of the maximized log-likelihood.}
\item{nPar}{Number of estimated parameters.}
\item{N}{Number of observations.}
\item{V}{Number of variables, corresponding to the number of nodes in the graph.}
\item{iter}{Number of iterations for the algorithm to converge.}
}

\details{
The function estimates a Gaussian graphical model given the graph association structure provided in input by the binary adjacency matrix. In the adjacency matrix, a zero entry corresponds to two variables being independent, marginally or conditionally according to the model.

If \code{model = "covariance"}, a \strong{Gaussian covariance graph model} is estimated, and the joint distribution of the \eqn{V}{V} dimensional vector of variables \eqn{X}{X} is parameterized in terms of the covariance matrix \eqn{\Sigma}{Sigma}. It is assumed:
\deqn{X \sim \mathcal{N}(\mu, \Sigma) \quad \Sigma \in C^{+}_G(A)}{X ~ N(mu, Sigma) with Sigma in C_G(A) }
where \eqn{C^{+}_G(A)}{C_G(A)} is the collection of sparse positive definite matrices whose zero patterns are given by graph \eqn{G}{G} represented by the adjacency matrix \eqn{A}{A}. In this type of model, the graph/adjacency matrix corresponds to \emph{marginal independence} constraints among the variables, i.e. the variables associated to two non-connected edges in the graph are marginally independent. As a result, the covariance matrix \code{sigma} is estimated to be sparse according to the graph.

If \code{model = "concentration"}, a \strong{Gaussian concentration graph model} is estimated, and the joint distribution of the \eqn{V}{V} dimensional vector of variables \eqn{X}{X} is parameterized in terms of the concentration matrix (inverse covariance or precision matrix) \eqn{\Omega}{Omega}. It is assumed:
\deqn{X \sim \mathcal{N}(\mu, \Omega) \quad \Omega \in C^{+}_G(A)}{X ~ N(mu, Omega) with Omega in C_G(A) }
where \eqn{C^{+}_G(A)}{C_G(A)} is the collection of sparse positive definite matrices whose zero patterns are given by graph \eqn{G}{G} embedded in the adjacency matrix \eqn{A}{A}. For this type of model, the graph/adjacency matrix corresponds to \emph{conditional independence} constraints among the variables, i.e. the variables associated to two non-adjacent edges in the graph are conditionally independent given their common neighbors. It results in the concentration matrix \code{omega} being estimated to be sparse according to the structure of the graph.

Note that conditional independence does not imply marginal independence, and marginal independence does not imply conditional independence, therefore a sparse concentration matrix and a sparse covariance matrix do not necessarily match; See Whittaker (1990).

The Gaussian covariance graph model is estimated using the \emph{Iterative Conditional Fitting} algorithm by Chaudhuri et al. (2007), while the Gaussian concentration graph model is estimated using the algorithm by Hastie et al. (2009).

Bayesian regularization is performed by means of a conjugate prior on the covariance/concentration matrix, similarly to what described in Fop et al. (2018). In the case of covariance graph model, an Inverse-Wishart distribution is used as a prior for \eqn{\Sigma}{Sigma}, while a Wishart distribution is used for \eqn{\Omega}{Omega} in the case of a concentration graph model. Regularization can be useful when the number of variables is larger than the number of observations.
}



\examples{

# Gaussian covariance graph model
data(mtcars)
x <- mtcars[,c(1,3:7)]
R <- cor(x)
#
# model where variables with correlation less than 0.5 are marginally independent
graph <- ( abs(R) < 0.5 )*1
diag(graph) <- 0
fit1 <- fitGGM(data = x, graph = graph)
plot(fit1)


# Gaussian concentration graph model
data(swiss)
#
# fit a conditional independence model:
V <- ncol(swiss)
graph <- matrix( c(0,1,0,1,1,1,
                   1,0,1,1,0,0,
                   0,1,0,1,1,0,
                   1,1,1,0,1,0,
                   1,0,1,1,0,0,
                   1,0,0,0,0,0), V,V, byrow = TRUE )
fit2 <- fitGGM(swiss, graph = graph, model = "concentration")
plot(fit2)


\dontrun{

data(marks, package = "ggm")
#
# the conditional independence model of Whittaker (1990), pag. 6
V <- ncol(marks)
graph <- matrix( c(0,1,1,0,0,
                   1,0,1,0,0,
                   1,1,0,1,1,
                   0,0,1,0,1,
                   0,0,1,1,0), V,V, byrow = TRUE )
fit3 <- fitGGM(marks, graph = graph, model = "concentration")
plot(fit3)
}

}

\references{
Chaudhuri, S., Drton M., and Richardson, T. S. (2007). Estimation of a covariance matrix with zeros. \emph{Biometrika}, 94(1), 199-216.

Fop, M., Murphy, T.B., and Scrucca, L. (2018). Model-based clustering with sparse covariance matrices. \emph{Statistics and Computing}. To appear.

Hastie, T., Tibshirani, R. and Friedman, J. (2009). \emph{The Elements of Statistical Learning}. Springer.

Whittaker, J. (1990). \emph{Graphical Models in Applied Multivariate Statistics}. Wiley.
}

\seealso{
\code{\link{plot.fitGGM}}
}
