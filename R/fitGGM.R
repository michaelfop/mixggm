#
#======================== Fit Gaussian graphical model
#
#

fitGGM <- function(data = NULL,
                   S = NULL, N = NULL,
                   graph,
                   model = c("covariance", "concentration"),    # inverse covariance to be implemented - omega # LSTODO: ??
                   start = NULL,
                   ctrlICF = controlICF(),
                   regularize = FALSE,
                   regHyperPar = NULL,
                   verbose = FALSE, ...)
{
  call <- match.call()
  if ( all(is.null(data), is.null(S)) ) 
    stop("We need some data to estimate a model! Please input 'data' or 'S' and 'N")
  if ( is.null(S) & !is.null(data) ) 
  {
    data <- data.matrix(data)
    N <- nrow(data)
    S <- cov(data)*(N-1)/N
  } else 
  if ( is.null(N) & is.null(data) ) 
    stop("You need to provide the sample size 'N' in input if don't supply 'data'")

  if(missing(graph))
    stop("'graph' argument is missing. Please provide a square symmetric binary adjacency matrix corresponding to the association structure of the graph.")
  graph <- as.matrix(graph)
  if ( !isSymmetric(graph) ) 
    stop ("'graph' must be a square symmetric binary adjacency matrix")
  # if ( any(diag(graph) != 0) ) 
  #   stop ("'graph' must be a square symmetric binary adjacency matrix with null diagonal")
  # LSTODO: forced to be binary with 0s along the diagonal
  diag(graph) <- 0
  graph[abs(graph) > 0] <- 1
  
  model <- match.arg(model, choices = eval(formals(fitGGM)$model))

  V <- ncol(S)
  if ( is.null(colnames(S) ) ) colnames(S) <- paste0("V", 1:V)
  varnames <- colnames(S)
  S <- as.matrix(S)
  nPar <- sum(graph)/2
  tot <- choose(V, 2)

  #### TODO: regularization for inverse covariance
  if ( regularize ) {
    if ( is.null(regHyperPar) ) {
      psi <- V + 2
      S <- if ( V > N ) {
        ( diag(diag(S)) + S*N ) / (psi + N + V + 1)
      } else ( S + S*N ) / (psi + N + V + 1)
    } else {
      psi <- regHyperPar$psi
      if ( !inherits(regHyperPar, "EM") ) S <- ( regHyperPar$scale + S*N ) / (regHyperPar$psi + N + V + 1)
      # if the function is not used in 'mixGGM' we compute the regularized S
      # if the function is used in 'mixGGM', regularized S is provided in input
    }
  } else {
    psi <- scale <- 0
  }

  if( min( eigen(S, only.values = TRUE)$val ) < sqrt(.Machine$double.eps)/10 )
    stop("Covariance matrix is not positive definite")

  # the graph is complete ....................................................
  if ( nPar == tot ) {
    sigma <- S
    if ( regularize ) N <- N + psi + V + 1
    lk <- profileloglik(sigma, S, N)
    dimnames(sigma) <- dimnames(lk$omega) <- dimnames(graph) <- list(varnames, varnames)
    res <- list(sigma = sigma, omega = lk$omega, graph = graph, model = model,
                loglik = lk$loglik, nPar = nPar + V, V = V, iter = 1)
    class(res) <- "fitGGM"
    return(res)
  }
  # the graph is fully disconnected ..........................................
  if ( nPar == 0 ) {
    sigma <- diag( diag(S) )
    dimnames(sigma) <- list(varnames, varnames)
    if ( regularize ) N <- N + psi + V + 1
    lk <- profileloglik(sigma, S, N)
    dimnames(lk$omega) <- dimnames(graph) <- list(varnames, varnames)
    res <- list(sigma = sigma, omega = lk$omega, graph = graph, model = model,
                loglik = lk$loglik, nPar = nPar + V, V = V, iter = 1)
    class(res) <- "fitGGM"
    return(res)
  }

  # get spouses and non spouses .............................................
  # SP <- NS <- SP2 <- list()
  # for ( i in 1:V ) {
  #   SP[[i]] <- which(graph[,i] == 1) - 1
  #   SP2[[i]] <- ifelse( SP[[i]] > i-1, SP[[i]] - 1, SP[[i]] )
  #   NS[[i]] <- setdiff(which(graph[,i] == 0), i) - 1
  # }
  # numSpo <- sapply(SP, length) != 0
  # nonTrivial <- which( numSpo != 0 )
  # noSpo <- which(numSpo == 0)
  # spouse <- findspouse(graph)

  # initialization ...........................................................
  if ( is.null(start) ) {    #### only used for covariance model
    sigma <- diag( diag(S) )
    # sigma <- S      # better?
  } else {
    temp <- diag(start)
    start[graph == 0] <- 0
    diag(start) <- temp
    sigma <- as.matrix(start)
    if ( min( eigen(sigma, only.values = TRUE)$values ) <= 0 ) 
      sigma <- diag( diag(S) )
  }
  #...........................................................................

  # icf ......................................................................
  out <- switch(model,
                covariance = icf(sigma, S, graph, N, 
                                 ctrlICF$tol, ctrlICF$itMax, 
                                 verbose, regularize, psi),
                concentration = conggm(S, graph, N, 
                                       ctrlICF$tol, ctrlICF$itMax, 
                                       verbose)
                )

  dimnames(out$sigma) <- dimnames(out$omega) <- list(varnames, varnames)
  res <- list(call = call, 
              model = model, graph = graph, N = N, V = V,
              loglik = out$loglik, iter = out$it, nPar = nPar + V, 
              sigma = out$sigma, omega = out$omega)
  class(res) <- "fitGGM"
  return(res)
}

print.fitGGM <- function(x, ...)
{
  if(!is.null(cl <- x$call))
  { 
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\n'fitGGM' object containing:","\n")
  print(names(x)[-1])
  invisible()
}

# print.fitGGM <- function(x, ...)
# {
#   cat("\n")
#   txt <- paste(" ", "Gaussian", x$model, "graph model", "\n")
#   cat(txt)
#   cat("    for", ifelse(x$model == "covariance", "marginal", "conditional"), "independence", "\n")
#   sep <- paste0(rep("=", max(nchar(txt)) + 1),
#                 collapse = "")
#   cat(sep, "\n")
#   cat( paste0("  ", "N. dependence parameters: ", x$nPar - x$V) )
#   cat("\n")
#   cat( paste0("  ", "Log-likelihood: ", round(x$loglik, 2), "\n") )
#   if ( !is.null(x$penalty) ) {
#     cat( paste0("  Penalized log-likelihood: ", round(x$loglikPen, 2), "\n") )
#     cat( paste0("  Penalty: ", x$penalty, "\n") )
#     cat( paste0("  Search: ", x$search, "\n") )
#   }
# }

