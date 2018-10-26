#
#======================== Graph structure GA search function
#
# Functions for performing graph search using genetic algorithms (poackage GA)
# in mixtures of Gaussian graphical models.

fitnessGGM <- function(string, par)
  # fitness function to be passed in the genetic algorithm
{
  N <- par$N

  # adjacency matrix ------------------------------------------------

  ########## THIS IS WRONG
  # adjMat <- matrix( 0, par$V, par$V )
  # adjMat[upper.tri(adjMat)] <- string
  # adjMat <- adjMat + t(adjMat)
  # dimnames(adjMat) <- list(colnames(par$S), colnames(par$S))
  #
  adjMat <- matrix( 0, par$V, par$V )
  adjMat[lower.tri(adjMat)] <- string
  adjMat <- adjMat + t(adjMat)
  dimnames(adjMat) <- list(colnames(par$S), colnames(par$S))

  # starting matrix for icf algorithm
  # start <- par$S/log2(par$V)
  start <- par$S
  start[ adjMat == 0 ] <- 0
  diag(start) <- diag(par$S)
  #------------------------------------------------------------------

  # fit covariance graph --------------------------------------------
  fit <- fitGGM(data = NULL, graph = adjMat, S = par$S, N = par$Nk, model = par$model, start = start,
                ctrlIcf = par$ctrlIcf, regularize = par$regularize, regHyperPar = par$regHyperPar)
  #------------------------------------------------------------------

  # compute criterion -----------------------------------------------
  crit <- fit$loglik - par$penalty(adjMat, beta = par$beta)
  #------------------------------------------------------------------

  return( crit )
}



searchGGMGA <- function(S, N, model = c("covariance", "concentration"),
                        pro = NULL, start = NULL,
                        regularize = FALSE, regHyperPar = NULL,
                        penalty = graphPenalty(), beta = NULL,
                        ctrlGa = ctrlGA(), ctrlIcf = ctrlICF(),
                        parallel = FALSE, verbose = FALSE)
  # Function for searching the best structure of a graph using genetic algorithms.
{
  V <- ncol(S)
  # if ( is.null(penPar) ) penPar <- list(N = N)
  model <- match.arg( model, c("covariance", "concentration") )

  # we need the weighted N in the profile loglikelihood for the EM algorithm
  if ( is.null(pro) ) Nk <- N else Nk <- N*pro

  if ( is.null(dimnames(S)) ) {
    varnames <- paste0("V", 1:V)
    dimnames(S) <- list(varnames, varnames)
  } else varnames <- colnames(S)

  allPairs <- combn( 1:V, 2, simplify = FALSE )
  nPairs <- length(allPairs)
  allCn <- sapply( allPairs, function(p) paste(p, collapse = "-") )

  # start parallel computations-----------------------------------------------
  if ( !inherits(parallel, "cluster") ) {    # if parallel not already started then
    if ( parallel | is.numeric(parallel) ) {
      parallel <- GA::startParallel(parallel)
      # cluster will be closed in 'ga' function
    }
  } else {  # if parallel already started
    parallel <- attr(parallel, "cluster")
    # else cluster will be closed in 'emCovGraph' function
  }
  #---------------------------------------------------------------------------

  # initialization -----------------------------------------------------------
  if ( is.null(start) ) {
    # provide subsets on the corner of the search space and on the correlation matrix
    R <- if ( model == "covariance" ) cov2cor(S) else cov2cor( solve(S) )
    corr <- abs( R[lower.tri(R)] )
    suggestions <- rbind( rep(0, nPairs), rep(1, nPairs),
                          ifelse( corr > median(corr), 1, 0 ),
                          ifelse( corr > 0.25, 1, 0 ), ifelse( corr > 0.3, 1, 0 ), ifelse( corr > 0.35, 1, 0 ),
                          ifelse( corr > 0.4, 1, 0 ), ifelse( corr > 0.45, 1, 0 ), ifelse( corr > 0.5, 1, 0 ),
                          ifelse( corr > 0.55, 1, 0 ), ifelse( corr > 0.6, 1, 0 ), ifelse( corr > 0.65, 1, 0 ),
                          ifelse( corr > 0.7, 1, 0 ), ifelse( corr > 0.75, 1, 0 ), ifelse( corr > 0.8, 1, 0 ),
                          ifelse( corr > 0.85, 1, 0 ), ifelse( corr > 0.9, 1, 0 ), ifelse( corr > 0.95, 1, 0 ))
  } else {
    zeroOne <- identical( range(start), c(0,1) ) | identical( range(start), c(0,0) )
    if ( !is.matrix(start) | !zeroOne )
      stop("You have to provide a proper adjacency matrix")
    suggestions <- matrix( start[lower.tri(start)], nrow = 1 )
  }
  #-------------------------------------------------------------------------

  # search the best graph --------------------------------------------------
  # parameters list
  par <- list(N = N, Nk = Nk, V = V, S = S, pro = pro, allPairs = allPairs,
              penalty = penalty, beta = beta, model = model,
              regularize = regularize, regHyperPar = regHyperPar, ctrlIcf = ctrlIcf)

  # memoise fitness function
  mfitness <- memoise::memoise(fitnessGGM)

  GA <- GA::ga(type = "binary",
               fitness = mfitness, par = par,
               nBits = nPairs,
               names = allCn,
               suggestions = suggestions,
               popSize = ctrlGa$popSize, pcrossover = ctrlGa$pcrossover, pmutation = ctrlGa$pmutation,
               maxiter = ctrlGa$maxiter, run = ctrlGa$run, elitism = ctrlGa$elitism,
               parallel = parallel, monitor = verbose
  )

  # reset cached memoised fitness function
  memoise::forget(mfitness)
  #-------------------------------------------------------------------------

  # estimated graph and association matrix-----------------------------------
  # graph <- matrix(0, V, V)
  # graphPairs <- allPairs[ GA@solution[1,] == 1 ]
  # for ( set in graphPairs ) graph[ set[1], set[2] ] <- graph[ set[2], set[1] ] <- 1
  #
  # faster
  graph <- matrix( 0, V, V )
  graph[lower.tri(graph)] <- GA@solution[1,]
  graph <- graph + t(graph)
  dimnames(graph) <- list(colnames(par$S), colnames(par$S))

  start <- S
  start[ graph == 0 ] <- 0
  diag(start) <- diag(S)

  dimnames(graph) <- list(varnames, varnames)
  fit <- fitGGM(data = NULL, graph = graph, S = S, N = Nk, model = model, regularize = regularize,
                regHyperPar = regHyperPar, ctrlIcf = ctrlIcf, start = start)
  #------------------------------------------------------------------------

  out <- list(sigma = fit$sigma,
              omega = fit$omega,
              graph = graph,
              loglik = fit$loglik,
              nPar = fit$nPar,
              crit = GA@fitnessValue,
              penalty = attr(penalty, "type"),
              GA = GA)
  return(out)
}
