#
#======================== Mixture of Gaussian graphical models
#
#

mixGGM <- function( data, K = 1:3,
                    model = c("covariance", "concentration"),
                    search = c("step-forw", "step-back","ga"),
                    penalty = c("bic", "ebic", "erdos", "power"),
                    beta = NULL,
                    regularize = FALSE, regHyperPar = NULL,
                    ctrlEm = ctrlEM(), ctrlStep = ctrlSTEP(), ctrlGa = ctrlGA(), ctrlIcf = ctrlICF(),
                    keepAll = FALSE,
                    parallel = FALSE,
                    verbose = TRUE )
  # Automatic model selection usin BIC
{
  varnames <- colnames(data)
  N <- nrow(data)
  V <- ncol(data)
  data <- data.matrix(data)
  model <- match.arg( model, c("covariance", "concentration") )
  search <- match.arg( search,  c("step-forw", "step-back","ga") )
  # penalty <- match.arg( penalty, c("bic", "ebic", "erdos", "power") )
  penalty <- graphPenalty(penalty)
  if ( is.null(varnames) ) {
    varnames <- paste0("V", 1:V)
    colnames(data) <- varnames
  }

  # start parallel computations--------------------------------------
  if ( parallel | is.numeric(parallel) ) {
    parallel <- GA::startParallel(parallel)
    class(parallel) <- "cluster"
    # inherits(ctrlGa$parallel, "cluster")
  }
  on.exit( if ( parallel ) parallel::stopCluster(attr(parallel, "cluster")) )
  #-------------------------------------------------------------------------

  # initialize EM algorithm
  if ( !is.null(ctrlEm$subset) ) {
    hcInit <- mclust::hc(data[ctrlEm$subset,], modelName = "VVV", use = "VARS")
  } else hcInit <- mclust::hc(data, modelName = "VVV", use = "VARS")

  res <- list()
  BIC <- rep( NA, length(K) )

  if ( verbose ) {
    pbar <- txtProgressBar(min = 0, max = length(K), style = 3)
    on.exit( close(pbar) )
  }

  if ( verbose ) setTxtProgressBar(pbar, 0)
  for ( k in K ) {
    i <- match(k, K)
    temp <- try( emMixGGM(data, K = k, model = model, search = search,
                          penalty = penalty, beta = beta,
                          regularize = regularize, regHyperPar = regHyperPar,
                          ctrlGa = ctrlGa, ctrlStep = ctrlStep,
                          ctrlEm = ctrlEm, ctrlIcf = ctrlIcf,
                          hcInit = hcInit, parallel = parallel),
                 silent = !ctrlEm$printMsg )
    if ( class(temp) == "try.error" | is.na(temp$loglik) ) {
      res[[i]] <- temp
      BIC[i] <- NA
    } else {
      res[[i]] <- temp
      # if ( regularize ) temp$loglik <- temp$llk
      BIC[i] <- temp$loglikReg - 0.5*log(N)*temp$nPar[2]
      BIC[i] <- BIC[i]*2    # to be consistent with mclust, otherwise comment

      # # icl
      # cl <- map(temp$z)
      # ICL[i] <- BIC[i] + 2*sum( cl * ifelse(temp$z > 0, log(temp$z), 0) )
      # # ebic
      # gamma <- 1   # fixed for now
      # TOT <- choose(M, 2)
      # Jk <- apply(temp$graph, 3, sum)/2
      # SJk <- choose(TOT, Jk)
      # EBIC[i] <- BIC[i] - 2*gamma*log( sum(SJk) )
      #
    }
    if ( verbose ) setTxtProgressBar(pbar, i)
  }

  # best model
  if ( all(is.na(BIC)) ) {
    h <- min(K)
    out <- temp[[h]]
    out$BIC <- out$bic <- NA
    out$data <- data
  } else {
    best <- which.max( BIC )
    fit <- res[[ best ]]
    out <- fit
    out$BIC <- BIC
    out$bic <- BIC[best]
    out$data <- data
    out$search <- search
    out$model<- model
  }

  # out <- out[c(1,2,16,5,3,4,6:9,13,12,10,11,15,14)]
  out <- out[ c(1:11, 14,13,15, 17,12,16)]
  out$keepAll <- if ( keepAll ) res else NULL
  out <- structure( out, control = list(EM = ctrlEm, STEP = ctrlStep, GA = ctrlGa, ICF = ctrlIcf) )
  class(out) <- "mixGGM"
  return(out)
}


print.mixGGM <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Mixture of Gaussian", x$model, "graph models", "\n")
  cat(txt)
  sep <- paste0(rep("=", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat( paste0("  ", "K = ", x$K, "\n") )
  cat( paste0("  ", "N. dependence parameters: ", x$nPar[1], "\n") )
  cat( paste0("  ", "Log-likelihood: ", round(x$loglik, 2), "\n") )
  cat( paste0("  Penalized log-likelihood: ", round(x$loglikPen, 2), "\n") )
  cat( paste0("  Penalty: ", x$penalty, "\n") )
  cat( paste0("  Search: ", x$search, "\n") )
}

