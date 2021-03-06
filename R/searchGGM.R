#
#======================== Graph structure search and model estimation
#
#

searchGGM <- function(data = NULL,
                      S = NULL, n = NULL,
                      model = c("covariance", "concentration"),
                      search = c("step-forw", "step-back", "ga"),   # lasso coming soon
                      penalty = c("bic", "ebic", "erdos", "power"),
                      beta = NULL,
                      start = NULL,
                      regularize = FALSE, regHyperPar = NULL,
                      ctrlSTEP = controlSTEP(), ctrlGA = controlGA(), ctrlICF = controlICF(),
                      parallel = FALSE,
                      verbose = FALSE, ...)
{
  if ( all(is.null(data), is.null(S)) ) stop("We need some data to estimate a model! Please input 'data' or 'S' and 'n")
  if ( is.null(S) & !is.null(data) ) {
    n <- nrow(data)
    S <- cov(data)*(n-1)/n
  } else if ( is.null(n) & is.null(data) ) stop("You need to provide the sample size 'n' in input if don't supply 'data'")

  V <- ncol(S)
  model <- match.arg(model, c("covariance", "concentration"))
  search <- match.arg(search, c("step-forw", "step-back", "ga"))
  penalty <- graphPenalty(penalty)

  forEM <- list(...)            # grab stuff for EM algorithm
  if ( length(forEM) < 1 ) {
    forEM <- list(pro = NULL, occam = NULL)
  }

  out <- switch(search,
                "step-back" = searchGGMStepwise_b(S = S, n = n, model = model,
                                                  penalty = penalty, beta = beta, start = start,
                                                  regularize = regularize, regHyperPar = regHyperPar,
                                                  ctrlSTEP = ctrlSTEP, ctrlICF = ctrlICF,
                                                  parallel = parallel, verbose = verbose,
                                                  pro = forEM$pro, occam = forEM$occam),
                "step-forw" = searchGGMStepwise_f(S = S, n = n, model = model,
                                                  penalty = penalty, beta = beta, start = start,
                                                  regularize = regularize, regHyperPar = regHyperPar,
                                                  ctrlSTEP = ctrlSTEP, ctrlICF = ctrlICF,
                                                  parallel = parallel, verbose = verbose,
                                                  pro = forEM$pro, occam = forEM$occam),
                "ga" = searchGGMGA(S = S, n = n, model = model,
                                   penalty = penalty, beta = beta, start = start,
                                   regularize = regularize, regHyperPar = regHyperPar,
                                   ctrlGA = ctrlGA, ctrlICF = ctrlICF,
                                   parallel = parallel, verbose = verbose,
                                   pro = forEM$pro)
  )
  res <- list(sigma = out$sigma, omega = out$omega, graph = out$graph, model = model, loglikPen = out$crit,
              loglik = out$loglik, nPar = out$nPar, n = n, V = V, penalty = attr(penalty, "type"), search = search)
  if ( search == "ga" ) res$GA <- out$GA
  class(res) <- "fitGGM"
  return(res)
}
