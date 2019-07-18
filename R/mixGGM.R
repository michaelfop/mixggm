#
#======================== Mixture of Gaussian graphical models
#
#

mixGGM <- function(data, K = 1:3,
                   model = c("covariance", "concentration"),
                   search = c("step-forw", "step-back","ga"),
                   penalty = c("bic", "ebic", "erdos", "power"),
                   beta = NULL,
                   regularize = FALSE, 
                   regHyperPar = NULL,
                   ctrlEM = controlEM(), 
                   ctrlSTEP = controlSTEP(), 
                   ctrlGA = controlGA(), 
                   ctrlICF = controlICF(),
                   keepAll = FALSE,
                   parallel = FALSE,
                   verbose = interactive() )
{
  call <- match.call()
  data <- data.matrix(data)
  model   <- match.arg(model, choices = eval(formals(mixGGM)$model))
  search  <- match.arg(search, choices = eval(formals(mixGGM)$search))
  penalty <- match.arg(penalty, choices = eval(formals(mixGGM)$penalty))
  penalty <- graphPenalty(penalty)
  
  varnames <- colnames(data)
  if(is.null(varnames)) 
  {
    varnames <- paste0("V", 1:V)
    colnames(data) <- varnames
  }
  N <- nrow(data)
  V <- ncol(data)

  # start parallel computations--------------------------------------
  if ( parallel | is.numeric(parallel) ) {
    parallel <- GA::startParallel(parallel)
    class(parallel) <- "cluster"
    # inherits(ctrlGA$parallel, "cluster")
  }
  on.exit( if ( parallel ) parallel::stopCluster(attr(parallel, "cluster")) )
  #-------------------------------------------------------------------------

  # initialize EM algorithm
  if ( !is.null(ctrlEM$subset) ) {
    hcInit <- mclust::hc(data[ctrlEM$subset,], modelName = "VVV", use = "VARS")
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
                          ctrlGA = ctrlGA, ctrlSTEP = ctrlSTEP,
                          ctrlEM = ctrlEM, ctrlICF = ctrlICF,
                          hcInit = hcInit, parallel = parallel),
                 silent = !ctrlEM$printMsg )
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
  # out <- out[ c(1:11, 14,13,15, 17,12,16)]
  out <- out[c(15,17,16,12,13,3,4,6,7,8,14,5,9,1,2,10,11)]
  out$keepAll <- if ( keepAll ) res else NULL

  out <- c(call = call, out)
  out$control <- list(EM = ctrlEM, 
                      STEP = ctrlSTEP, 
                      GA = ctrlGA, 
                      ICF = ctrlICF)
  class(out) <- "mixGGM"
  return(out)
}

print.mixGGM <- function(x, ...)
{
  if(!is.null(cl <- x$call))
  { 
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\n'mixGGM' object containing:","\n")
  print(names(x)[-1])
  invisible()
}

# print.mixGGM <- function(x, ...)
# {
#   cat("\n")
#   txt <- paste(" ", "Mixture of Gaussian", x$model, "graph models", "\n")
#   cat(txt)
#   sep <- paste0(rep("=", max(nchar(txt)) + 1),
#                 collapse = "")
#   cat(sep, "\n")
#   cat( paste0("  ", "K = ", x$K, "\n") )
#   cat( paste0("  ", "N. dependence parameters: ", x$nPar[1], "\n") )
#   cat( paste0("  ", "Log-likelihood: ", round(x$loglik, 2), "\n") )
#   cat( paste0("  Penalized log-likelihood: ", round(x$loglikPen, 2), "\n") )
#   cat( paste0("  Penalty: ", x$penalty, "\n") )
#   cat( paste0("  Search: ", x$search, "\n") )
# }

summary.mixGGM <- function(object, 
                           graphs = TRUE, 
                           clusters = FALSE,
                           parameters = FALSE,
                           ...)
{
  out <- object[c("model", "penalty", "search", "N", "V", 
                  "loglik", "loglikPen", "loglikReg",
                  "K", "nPar", "bic", "graph", "parameters")]
  out$regularize <- (object$loglik != object$loglikReg)
  out$printGraphs <- graphs
  out$printClusters <- clusters
  out$printParameters <- parameters
  out$tabClassification <- table(factor(object$classification))
  out$search <- switch(out$search,
                       "step-forw" = "forward-stepwise",
                       "step-back" = "backward-stepwise",
                       "ga" = "genetic algorithm")
  out$penalty <- if(is.character(out$penalty)) out$penalty else "user defined"
  class(out) <- "summary.mixGGM"
  return(out)
}

print.summary.mixGGM <- function(x, digits = getOption("digits"), ...)
{
  cat(cli::rule(left = crayon::bold("Mixture of Gaussian graphical models"), 
                width = min(getOption("width"),60)), "\n\n")
  #
  cat(paste("Data dimensions =", x$N, "x", x$V, "\n"))
  cat(paste("Model           =", x$model, "\n"))
  cat(paste("Search          =", x$search, "\n"))
  cat(paste("Penalty         =", x$penalty, "\n\n"))
  #
  tab <- data.frame("l" = if(x$regularize) x$loglikReg else x$loglik,
                    "Pen. log-likelihood" = x$loglikPen, 
                    "df" = x$nPar[2], "BIC" = x$bic, 
                    row.names = "", check.names = FALSE)
  colnames(tab)[1] <- if(x$regularize) "Reg. log-likelihood" else "Log-likelihood"
  print(tab, digits = digits)
  #
  if(x$printGraphs)
  {
    cat("\nGraphs:\n") 
    for(k in 1:x$K)
    { 
      cat("[,,", k, "]", sep = "")
      g <- as.data.frame(x$graph[,,k])
      colnames(g) <- rep("", ncol(g))
      # print(ifelse(g == 0, ".", "*"), quote = FALSE)
      print(ifelse(g == 0, cli::symbol$circle, cli::symbol$circle_filled), quote = FALSE)
      # 
    }
  }
  #
  if(x$printClusters)
  { 
    cat("\nClustering table:")
    print(x$tabClassification, digits = digits)
  }
  #
  if(x$printParameters)
  { 
    cat("\nMixing probabilities:\n")
    print(x$parameters$tau, digits = digits)
    cat("\nMeans:\n")
    print(x$parameters$mu, digits = digits)
    if(x$model == "covariance") 
    { 
      cat("\nCovariances:\n") 
      for(k in 1:x$K)
      { 
        cat("[,,", k, "]\n", sep = "")
        print(x$parameters$sigma[,,k], digits = digits)
      }
    } else 
    {
      cat("\nConcentrations:\n")
      for(k in 1:x$K)
      { 
        cat("[,,", k, "]\n", sep = "")
        print(x$parameters$omega[,,k], digits = digits)
      }
    }
  }
  #
  invisible()
}
