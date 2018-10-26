#
#=============================== Utility control functions
#


ctrlICF <- function(tol = 1e-04, maxiter = 1e03)
  # icf parameters
{
  list( tol = tol, itMax = maxiter )
}


ctrlREG <- function(data, K,
                    scaleType = c("full", "fixed", "one", "diag"),
                    scale = NULL, psi = NULL)
  # hyperparameters for Bayesian regularization
{
  V <- ncol(data)
  N <- nrow(data)
  st <- match.arg( scaleType, c("full", "fixed", "one", "diag") )

  if ( is.null(psi) ) psi <- V + 2
  if ( is.null(scale) ) {
    VAR <- var(data)*(N-1)/N
    if ( N > V ) {
      scale <- switch( st,
                       full = VAR/( K^(2/V) ),
                       fixed = VAR/det(VAR)^(1/V) * (0.001/K)^(1/V),
                       one = VAR/det(VAR)^(1/V),
                       diag = diag( diag(VAR)/( K^(2/V) ) )
      )
    } else {
      scale <- diag( diag(VAR)/( K^(2/V) ) )
    }
  }
  out <- list(psi = psi, scale = scale)
  return(out)
}


ctrlEM <- function(tol = 1e-05, maxiter = 1e02, subset = NULL, printMsg = FALSE)
  # EM control parameters
{
  list(tol = tol, maxiter = maxiter, subset = subset, printMsg = printMsg)
}


ctrlGA <- function(popSize = 50, pcrossover = 0.8, pmutation = 0.1,
                   maxiter = 100, run = maxiter/2,
                   elitism = base::max(1, round(popSize*0.05)))
  # GA search parameters
{
  list( popSize = popSize, pcrossover = pcrossover, pmutation = pmutation,
        maxiter = maxiter, run = run, elitism = elitism)
}


ctrlSTEP <- function(occamAdd = Inf, occamRem = Inf, start = NULL)
  # stepwise search parameters
{
  list(occamAdd = occamAdd, occamRem = occamRem, start = start)
}


### no longer in use
# profileLogLik <- function(Sigma, S, N)
#   # Compute Gaussian profile log-likelihood
# {
#   V <- unique( c(dim(S), dim(Sigma)) )
#   if ( length(V) > 1 ) stop("Wrong input parameters!")
#   inv <- solve(Sigma)
#   val <- -N/2*determinant(Sigma)$modulus - N/2*sum( diag( crossprod(inv, S) ) )
#   return( list(val = val, inv = inv) )
# }
