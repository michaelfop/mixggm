#
#======================== Graph structure search and model estimation
#
#

searchGGM <- function(data = NULL,
                      S = NULL, N = NULL,
                      model = c("covariance", "concentration"),
                      search = c("step-forw", "step-back", "ga"),   # lasso coming soon
                      penalty = c("bic", "ebic", "erdos", "power"),
                      beta = NULL,
                      start = NULL,
                      regularize = FALSE, regHyperPar = NULL,
                      ctrlStep = ctrlSTEP(), ctrlGa = ctrlGA(), ctrlIcf = ctrlICF(),
                      parallel = FALSE,
                      verbose = FALSE, ...)
{
  if ( all(is.null(data), is.null(S)) ) stop("We need some data to estimate a model! Please input 'data' or 'S' and 'N")
  if ( is.null(S) & !is.null(data) ) {
    N <- nrow(data)
    S <- cov(data)*(N-1)/N
  } else if ( is.null(N) & is.null(data) ) stop("You need to provide the sample size 'N' in input if don't supply 'data'")

  V <- ncol(S)
  model <- match.arg(model, c("covariance", "concentration"))
  search <- match.arg(search, c("step-forw", "step-back", "ga"))
  penalty <- graphPenalty(penalty)

  forEM <- list(...)            # grab stuff for EM algorithm
  if ( length(forEM) < 1 ) {
    forEM <- list(pro = NULL, occam = NULL)
  }

  out <- switch(search,
                "step-back" = searchGGMStepwise_b(S = S, N = N, model = model,
                                                  penalty = penalty, beta = beta, start = start,
                                                  regularize = regularize, regHyperPar = regHyperPar,
                                                  ctrlStep = ctrlStep, ctrlIcf = ctrlIcf,
                                                  parallel = parallel, verbose = verbose,
                                                  pro = forEM$pro, occam = forEM$occam),
                "step-forw" = searchGGMStepwise_f(S = S, N = N, model = model,
                                                  penalty = penalty, beta = beta, start = start,
                                                  regularize = regularize, regHyperPar = regHyperPar,
                                                  ctrlStep = ctrlStep, ctrlIcf = ctrlIcf,
                                                  parallel = parallel, verbose = verbose,
                                                  pro = forEM$pro, occam = forEM$occam),
                "ga" = searchGGMGA(S = S, N = N, model = model,
                                   penalty = penalty, beta = beta, start = start,
                                   regularize = regularize, regHyperPar = regHyperPar,
                                   ctrlGa = ctrlGa, ctrlIcf = ctrlIcf,
                                   parallel = parallel, verbose = verbose,
                                   pro = forEM$pro)
  )
  res <- list(sigma = out$sigma, omega = out$omega, graph = out$graph, model = model, loglikPen = out$crit,
              loglik = out$loglik, nPar = out$nPar, N = N, V = V, penalty = attr(penalty, "type"), search = search)
  if ( search == "ga" ) res$GA <- out$GA
  class(res) <- "fitGGM"
  return(res)
}
