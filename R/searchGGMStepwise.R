#
#======================== Graph structure stepwise search function
#
# Functions for performing graph search using a forward stepwise algorithm
# in mixtures of Gaussian graphical models

startStepwise <- function(S, Nk, N, model, penalty, beta, regularize, regHyperPar, ctrlIcf)
  # initialize graph in stepwise search
{
  V <- ncol(S)
  nPairs <- choose(V, 2)
  # if ( is.null(penPar) ) penPar <- list(N = N)

  # provide subsets on the corner of the search space and on the correlation matrix
  R <- if ( model == "covariance" ) cov2cor(S) else cov2cor( solve(S) )
  corr <- abs( R[lower.tri(R)] )
  sugg <- rbind( rep(0, nPairs), # rep(1, nPairs),     # never suggest full model
                 ifelse( corr > mean(corr), 1, 0 ),
                 ifelse( corr > median(corr), 1, 0 ),  # half number of parameters
                 # ifelse( corr > 0.35, 1, 0 ),
                 ifelse( corr > 0.4, 1, 0 ), ifelse( corr > 0.45, 1, 0 ),
                 ifelse( corr > 0.5, 1, 0 ), ifelse( corr > 0.55, 1, 0 ), ifelse( corr > 0.6, 1, 0 ),
                 ifelse( corr > 0.65, 1, 0 ), ifelse( corr > 0.7, 1, 0 ), ifelse( corr > 0.75, 1, 0 ),
                 ifelse( corr > 0.8, 1, 0 ), ifelse( corr > 0.85, 1, 0 ), ifelse( corr > 0.9, 1, 0 ),
                 ifelse( corr > 0.95, 1, 0 ))
  sugg <- sugg[!duplicated(sugg),,drop = FALSE]
  ns <- nrow(sugg)
  A <- matrix(0, V,V)
  ut <- upper.tri(A)
  lt <- lower.tri(A)
  startIcf <- S/log2(V)    # starting matrix for icf algorithm
  crit <- rep(NA, ns)
  fit <- vector("list", ns)

  for ( s in 1:ns ) {
    A[lt] <- sugg[s,]
    A[ut] <- t(A)[ut]
    startIcf <- S/log2(V)
    startIcf[ A == 0 ] <- 0
    diag(startIcf) <- diag(S)
    fit[[s]] <- fitGGM(data = NULL, S = S, graph = A, model = model, N = Nk, start = startIcf, ctrlIcf = ctrlIcf,
                       regularize = regularize, regHyperPar = regHyperPar)
    crit[s] <- fit[[s]]$loglik - penalty(A, beta = beta)
  }

  ret <- which.max(crit)
  A[lt] <- sugg[ret,]
  A[ut] <- t(A)[ut]
  out <- list(adjMat = A, Sigma = fit[[ret]]$sigma, Omega = fit[[ret]]$omega, val = crit[ret])
  return(out)
}



searchGGMStepwise_f <- function(S, N, model = c("covariance", "concentration"),
                                pro = NULL, start = NULL,
                                regularize = FALSE, regHyperPar = NULL,
                                penalty = graphPenalty(), beta = NULL,
                                ctrlStep = ctrlSTEP(), ctrlIcf = ctrlICF(),
                                parallel = FALSE, verbose = FALSE, occam = NULL)
  # Function for searching the best structure of a graph using a forward stepwise search.
{
  k <- NULL # avoid 'no visible binding for global variable'

  V <- ncol(S)
  nPairs <- choose(V, 2)
  # if ( is.null(penPar) ) penPar <- list(N = N)

  # we need the weighted N in the profile loglikelihood
  Nk <- if ( is.null(pro) ) N else N*pro
  model <- match.arg( model, c("covariance", "concentration") )

  if ( is.null(dimnames(S)) ) {
    cnames <- paste0("V", 1:V, sep = " ")
    dimnames(S) <- list(cnames, cnames)
  } else cnames <- colnames(S)

  # start parallel computations-----------------------------------------------
  if ( !inherits(parallel, "cluster") ) {    # if parallel not already started then
    if ( parallel | is.numeric(parallel) ) {
      parallel <- GA::startParallel(parallel)
      # cluster will be closed here
      on.exit( if (parallel) parallel::stopCluster(attr(parallel, "cluster")) )
    }
  } else {  # if parallel already started
    parallel <- attr(parallel, "cluster")
    # else cluster will be closed in 'emCovGraph' function
  }
  `%DO%` <- if ( is.logical(parallel) ) {
    if ( parallel ) foreach::`%dopar%` else foreach::`%do%`
  } else foreach::`%dopar%`
  #---------------------------------------------------------------------------

  # initialization ---------------------------------------------------------
  if ( is.null(start) ) {
    startStep <- startStepwise(S, Nk, N, model, penalty, beta, regularize, regHyperPar, ctrlIcf)
    adjMat0 <- startStep$adjMat
    Sigma <- startStep$Sigma
    Omega <- startStep$Omega
    critVal <- startStep$val
  } else {
    if ( is.null(attr(start, "critOut")) ) attr(start, "critOut") <- NA
    if ( is.na(attr(start, "critOut")) ) {
      zeroOne <- identical( range(start), c(0,1) ) | identical( range(start), c(0,0) )
      if ( !is.matrix(start) | !zeroOne )
        stop("You have to provide a proper adjacency matrix")
      adjMat0 <- start
      fit <- fitGGM(graph = adjMat0, S = S, N = Nk, model = model, start = S, ctrlIcf = ctrlIcf,
                    regularize = regularize, regHyperPar = regHyperPar)
      Sigma <- fit$sigma
      Omega <- fit$omega
      npar <- sum(adjMat0)/2
      critVal <- fit$loglik - penalty(adjMat0, beta = beta)
    } else {
      #----------- to be used within the em algorithm (also if start is provided in the em)
      # Sigma <- attr(start, "sigma")
      critVal <- attr(start, "critOut")
      adjMat0 <- start
      attributes(adjMat0)[c("sigma", "critOut")] <- NULL
      fit <- fitGGM(data = NULL, graph = adjMat0, S = S, N = Nk, model = model, start = S, ctrlIcf = ctrlIcf,
                    regularize = regularize, regHyperPar = regHyperPar)
      Sigma <- fit$sigma
      Omega <- fit$omega  # update matrices
      # critVal <- fit$loglik + penalty(adjMat0, par = penPar)  # update critVal after z estimation
      # critVal <- max( critVal, fit$loglik + penalty(adjMat0, par = penPar) )
    }
  }
  # ------------------------------------------------------------------------

  # algorithm parameters and start -----------------------------------------
  adjMat <- adjMat0
  crit <- TRUE
  niter <- 0

  occamAdd <- ctrlStep$occamAdd
  occamRem <- ctrlStep$occamRem
  inWindowAdd <- inWindowRem <- TRUE
  upDiag <- upper.tri(adjMat)*adjMat
  TODROP <- if ( is.null(occam$TODROP) ) which(upDiag != 0, arr.ind = TRUE) else occam$TODROP    # edges to remove
  lowDiag <- lower.tri(adjMat, diag = TRUE)
  upDiag[lowDiag] <- -1
  TOADD <- if ( is.null(occam$TOADD) ) which(upDiag == 0, arr.ind = TRUE) else occam$TOADD     # edges to add
  nodein <- nodeout <- nPairs + 1
  added <- removed <- FALSE
  # ------------------------------------------------------------------------

  while ( crit ) {

    # add an edge --------------------------------------------------------
    toAdd <- if ( length(TOADD) > 0 ) TOADD[-nodeout,,drop = FALSE] else NULL          # don't add the last edge just removed
    nr <- if ( !is.null(toAdd) ) nrow(toAdd) else 0
    #
    if ( nr > 0 ) {
      # if nr !> 0 there are no edges to add
      outVal <- foreach::foreach( k = 1:nr, .multicombine = TRUE, .errorhandling = "remove" ) %DO% {
        i <- toAdd[k,1]
        j <- toAdd[k,2]
        a <- adjMat
        a[i,j] <- a[j,i] <- 1
        #
        fit <- fitGGM(data = NULL, graph = a, S = S, N = Nk, model = model, start = Sigma, ctrlIcf = ctrlIcf,
                      regularize = regularize, regHyperPar = regHyperPar)
        #
        val <- fit$loglik - penalty(a, beta = beta)
        out <- list( fit$sigma, fit$omega, val )
        # the first element of the list 'out' contains the estimated covariance
        # matrices, the second one contains the criterion values
        return(out)
      }
      valMat <- sapply(outVal, "[[", 3)
      best <- which.max( valMat )
      sel <- toAdd[best,, drop = FALSE]
      selRev <- rev(sel)
      attributes(selRev) <- attributes(sel)
      Max <- valMat[best]

      if ( Max > critVal ) {
        added <- TRUE
        TODROP <- rbind(TODROP, sel)
        nodein <- nrow(TODROP)
        TOADD <- TOADD[-best,,drop = FALSE]
        adjMat[sel] <- adjMat[selRev] <- 1
        critVal <- Max
        Sigma <- outVal[[best]][[1]]
        Omega <- outVal[[best]][[2]]
        if ( verbose )
          cat( " Added", paste(cnames[c(sel)], collapse = "--"),
               "  CRIT =", Max, "\n")
        #
        # edges to add within the Occam's window
        diffVal <- valMat[-best] - critVal
        inWindowAdd <- diffVal > -occamAdd
        TOADD <- if ( length(inWindowAdd) > 0 ) TOADD[inWindowAdd,,drop = FALSE] else inWindowAdd
      } else {
        if ( sum(adjMat) == 0 ) break     # the selected model is the diagonal one
        added <- FALSE
        nodein <- nPairs + 1
        #
        diffVal <- valMat - critVal
        inWindowAdd <- diffVal > -occamAdd
        TOADD <- if ( removed ) {
          rbind(toAdd[inWindowAdd,], TOADD[nodeout,,drop = FALSE])
        } else toAdd[inWindowAdd,,drop = FALSE]
      }

      if ( sum(adjMat0) == 0 & niter < 1 ) {
        # if we started from the diagonal model we repeat an adding step
        niter <- 1
        next
      }

    } else {    # if nr is not > 0
      added <- FALSE
      nodein <- nPairs + 1
    }
    # -----------------------------------------------------------

    # remove one edge -------------------------------------------
    toDrop <- if ( length(TODROP) > 0 ) TODROP[-nodein,,drop = FALSE] else NULL    # don't remove the last edge just added
    nr <- if ( !is.null(toDrop) ) nrow(toDrop) else 0
    #
    if ( nr > 0 & (added | removed | niter < 1) ) {           # if nr !> 0 there are no edges to remove
      # if two consecutive steps not accepted, stop
      outVal <- foreach::foreach( k = 1:nr, .multicombine = TRUE, .errorhandling = "remove" ) %DO% {
        i <- toDrop[k,1]
        j <- toDrop[k,2]
        a <- adjMat
        a[i,j] <- a[j,i] <- 0
        #
        fit <- fitGGM(data = NULL, graph = a, S = S, N = Nk, model = model, start = Sigma, ctrlIcf = ctrlIcf,
                      regularize = regularize, regHyperPar = regHyperPar)
        #
        val <- fit$loglik - penalty(a, beta = beta)
        out <- list( fit$sigma, fit$omega, val )
        return(out)
      }
      valMat <- sapply(outVal, "[[", 3)
      best <- which.max( valMat )
      sel <- toDrop[best,, drop = FALSE]
      selRev <- rev(sel)
      attributes(selRev) <- attributes(sel)
      Max <- valMat[best]

      if ( Max >= critVal ) {
        removed <- TRUE
        TOADD <- rbind(TOADD, sel)
        nodeout <- nrow(TOADD)
        TODROP <- TODROP[-best,,drop = FALSE]
        adjMat[sel] <- adjMat[selRev] <- 0
        critVal <- Max
        Sigma <- outVal[[best]][[1]]
        Omega <- outVal[[best]][[2]]
        if ( verbose )
          cat( " Removed", paste(cnames[c(sel)], collapse = "--"),
               "  CRIT =", Max, "\n")
        #
        # edges to remove within the Occam's window
        diffVal <- valMat[-best] - critVal
        inWindowRem <- diffVal > -occamRem
        TODROP <- if ( length(inWindowRem) > 0 ) TODROP[inWindowRem,,drop = FALSE] else inWindowRem
      } else {
        if ( sum(adjMat)/2 == nPairs ) break      # the selected model is the full one
        removed <- FALSE
        nodeout <- nPairs + 1
        #
        diffVal <- valMat - critVal
        inWindowRem <- diffVal > -occamRem
        TODROP <- if ( added ) {
          rbind(toDrop[inWindowRem,], TODROP[nodein,,drop = FALSE])
        } else toDrop[inWindowRem,,drop = FALSE]
      }

      #
    } else {    # if nr is not > 0
      removed <- FALSE
      nodeout <- nPairs + 1
    }
    # -----------------------------------------------------------

    # check -----------------------------------------------------
    niter <- niter + 1
    # print( niter )
    crit <- ( added == TRUE | removed == TRUE )
  } # while

  # if ( parallel ) stopCluster(cl)
  out <- list(sigma = Sigma,
              omega = Omega,
              graph = adjMat,
              loglik = critVal + penalty(adjMat, beta = beta),
              nPar = sum(adjMat)/2 + V,
              crit = critVal,
              penalty = attr(penalty, "type"),
              occam = list(TOADD = TOADD, TODROP = TODROP))
  return(out)
}


searchGGMStepwise_b <- function(S, N, model = c("covariance", "concentration"),
                                pro = NULL, start = NULL,
                                regularize = FALSE, regHyperPar = NULL,
                                penalty = graphPenalty(), beta = NULL,
                                ctrlStep = ctrlSTEP(), ctrlIcf = ctrlICF(),
                                parallel = FALSE, verbose = FALSE, occam = NULL)
  # Function for searching the best structure of a graph using a backward stepwise search.
{
  k <- NULL # avoid 'no visible binding for global variable'

  V <- ncol(S)
  nPairs <- choose(V, 2)
  # if ( is.null(penPar) ) penPar <- list(N = N)

  # we need the weighted N in the profile loglikelihood
  Nk <- if ( is.null(pro) ) N else N*pro
  model <- match.arg( model, c("covariance", "concentration") )

  if ( is.null(dimnames(S)) ) {
    cnames <- paste0("V", 1:V, sep = " ")
    dimnames(S) <- list(cnames, cnames)
  } else cnames <- colnames(S)

  # start parallel computations-----------------------------------------------
  if ( !inherits(parallel, "cluster") ) {    # if parallel not already started then
    if ( parallel | is.numeric(parallel) ) {
      parallel <- GA::startParallel(parallel)
      # cluster will be closed here
      on.exit( if (parallel) parallel::stopCluster(attr(parallel, "cluster")) )
    }
  } else {  # if parallel already started
    parallel <- attr(parallel, "cluster")
    # else cluster will be closed in 'emCovGraph' function
  }
  `%DO%` <- if ( is.logical(parallel) ) {
    if ( parallel ) foreach::`%dopar%` else foreach::`%do%`
  } else foreach::`%dopar%`
  #---------------------------------------------------------------------------

  # initialization ---------------------------------------------------------
  if ( is.null(start) ) {
    # start from fully connected graph
    adjMat0 <- matrix(1, V,V)
    diag(adjMat0) <- 0
    Sigma <- S
    val <- profileloglik(Sigma, S, N)
    Omega <- val$inv
    critVal <- val$loglik - penalty(adjMat0, beta = beta)
  } else {
    if ( is.null(attr(start, "critOut")) ) attr(start, "critOut") <- NA
    if ( is.na(attr(start, "critOut")) ) {
      zeroOne <- identical( range(start), c(0,1) ) | identical( range(start), c(0,0) )
      if ( !is.matrix(start) | !zeroOne )
        stop("You have to provide a proper adjacency matrix")
      adjMat0 <- start
      fit <- fitGGM(graph = adjMat0, S = S, N = Nk, model = model, start = S, ctrlIcf = ctrlIcf,
                    regularize = regularize, regHyperPar = regHyperPar)
      Sigma <- fit$sigma
      Omega <- fit$omega
      npar <- sum(adjMat0)/2
      critVal <- fit$loglik - penalty(adjMat0, beta = beta)
    } else {
      #----------- to be used within the em algorithm (also if start is provided in the em)
      # Sigma <- attr(start, "sigma")
      critVal <- attr(start, "critOut")
      adjMat0 <- start
      attributes(adjMat0)[c("sigma", "critOut")] <- NULL
      fit <- fitGGM(data = NULL, graph = adjMat0, S = S, N = Nk, model = model, start = S, ctrlIcf = ctrlIcf,
                    regularize = regularize, regHyperPar = regHyperPar)
      Sigma <- fit$sigma
      Omega <- fit$omega  # update matrices
      # critVal <- fit$loglik + penalty(adjMat0, par = penPar)  # update critVal after z estimation
      # critVal <- max( critVal, fit$loglik + penalty(adjMat0, par = penPar) )
    }
  }
  # ------------------------------------------------------------------------

  # algorithm parameters and start -----------------------------------------
  adjMat <- adjMat0
  crit <- TRUE
  niter <- 0

  occamAdd <- ctrlStep$occamAdd
  occamRem <- ctrlStep$occamRem
  inWindowAdd <- inWindowRem <- TRUE
  upDiag <- upper.tri(adjMat)*adjMat
  TODROP <- if ( is.null(occam$TODROP) ) which(upDiag != 0, arr.ind = TRUE) else occam$TODROP    # edges to remove
  lowDiag <- lower.tri(adjMat, diag = TRUE)
  upDiag[lowDiag] <- -1
  TOADD <- if ( is.null(occam$TOADD) ) which(upDiag == 0, arr.ind = TRUE) else occam$TOADD     # edges to add
  nodein <- nodeout <- nPairs + 1
  added <- removed <- FALSE
  # ------------------------------------------------------------------------

  while ( crit ) {

    # remove one edge -------------------------------------------
    toDrop <- if ( length(TODROP) > 0 ) TODROP[-nodein,,drop = FALSE]  else NULL # don't remove the last edge just added
    nr <- if ( !is.null(toDrop) ) nrow(toDrop) else 0
    #
    # if ( nr > 0 & (added | removed | niter < 1) ) {           # if nr !> 0 there are no edges to remove
    # if two consecutive steps not accepted, stop
    if ( nr > 0 ) {
      outVal <- foreach::foreach( k = 1:nr, .multicombine = TRUE, .errorhandling = "remove" ) %DO% {
        i <- toDrop[k,1]
        j <- toDrop[k,2]
        a <- adjMat
        a[i,j] <- a[j,i] <- 0
        #
        fit <- fitGGM(data = NULL, graph = a, S = S, N = Nk, model = model, start = Sigma, ctrlIcf = ctrlIcf,
                      regularize = regularize, regHyperPar = regHyperPar)
        #
        val <- fit$loglik - penalty(a, beta = beta)
        out <- list( fit$sigma, fit$omega, val )
        return(out)
      }
      valMat <- sapply(outVal, "[[", 3)
      best <- which.max( valMat )
      sel <- toDrop[best,, drop = FALSE]
      selRev <- rev(sel)
      attributes(selRev) <- attributes(sel)
      Max <- valMat[best]

      if ( Max >= critVal ) {
        removed <- TRUE
        TOADD <- rbind(TOADD, sel)
        nodeout <- nrow(TOADD)
        TODROP <- TODROP[-best,,drop = FALSE]
        adjMat[sel] <- adjMat[selRev] <- 0
        critVal <- Max
        Sigma <- outVal[[best]][[1]]
        Omega <- outVal[[best]][[2]]
        if ( verbose )
          cat( " Removed", paste(cnames[c(sel)], collapse = "--"),
               "  CRIT =", Max, "\n")
        #
        # edges to remove within the Occam's window
        diffVal <- valMat[-best] - critVal
        inWindowRem <- diffVal > -occamRem
        TODROP <- if ( length(inWindowRem) > 0 ) TODROP[inWindowRem,,drop = FALSE] else inWindowRem
      } else {
        if ( sum(adjMat)/2 == nPairs ) break      # the selected model is the full one
        removed <- FALSE
        nodeout <- nPairs + 1
        #
        diffVal <- valMat - critVal
        inWindowRem <- diffVal > -occamRem
        TODROP <- if ( added ) {
          rbind(toDrop[inWindowRem,], TODROP[nodein,,drop = FALSE])
        } else toDrop[inWindowRem,,drop = FALSE]
      }

      if ( niter < 1 ) {
        # we start with two adding steps
        niter <- 1
        next
      }

      #
    } else {    # if nr is not > 0
      removed <- FALSE
      nodeout <- nPairs + 1
    }
    # -----------------------------------------------------------

    # add an edge --------------------------------------------------------
    toAdd <- if ( length(TOADD) > 0 ) TOADD[-nodeout,,drop = FALSE] else NULL          # don't add the last edge just removed
    nr <- if ( !is.null(toAdd) ) nrow(toAdd) else 0
    #
    if ( nr > 0 & ( added | removed ) ) {
      # if nr !> 0 there are no edges to add
      # also if two consecutive steps not accepted, stop
      outVal <- foreach::foreach( k = 1:nr, .multicombine = TRUE, .errorhandling = "remove" ) %DO% {
        i <- toAdd[k,1]
        j <- toAdd[k,2]
        a <- adjMat
        a[i,j] <- a[j,i] <- 1
        #
        fit <- fitGGM(data = NULL, graph = a, S = S, N = Nk, model = model, start = Sigma, ctrlIcf = ctrlIcf,
                      regularize = regularize, regHyperPar = regHyperPar)
        #
        val <- fit$loglik - penalty(a, beta = beta)
        out <- list( fit$sigma, fit$omega, val )
        # the first element of the list 'out' contains the estimated covariance
        # matrices, the second one contains the criterion values
        return(out)
      }
      valMat <- sapply(outVal, "[[", 3)
      best <- which.max( valMat )
      sel <- toAdd[best,, drop = FALSE]
      selRev <- rev(sel)
      attributes(selRev) <- attributes(sel)
      Max <- valMat[best]

      if ( Max > critVal ) {
        added <- TRUE
        TODROP <- rbind(TODROP, sel)
        nodein <- nrow(TODROP)
        TOADD <- TOADD[-best,,drop = FALSE]
        adjMat[sel] <- adjMat[selRev] <- 1
        critVal <- Max
        Sigma <- outVal[[best]][[1]]
        Omega <- outVal[[best]][[2]]
        if ( verbose )
          cat( " Added", paste(cnames[c(sel)], collapse = "--"),
               "  CRIT =", Max, "\n")
        #
        # edges to add within the Occam's window
        diffVal <- valMat[-best] - critVal
        inWindowAdd <- diffVal > -occamAdd
        TOADD <- if ( length(inWindowAdd) > 0 ) TOADD[inWindowAdd,,drop = FALSE] else inWindowAdd
      } else {
        if ( sum(adjMat) == 0 ) break     # the selected model is the diagonal one
        added <- FALSE
        nodein <- nPairs + 1
        #
        diffVal <- valMat - critVal
        inWindowAdd <- diffVal > -occamAdd
        TOADD <- if ( removed ) {
          rbind(toAdd[inWindowAdd,], TOADD[nodeout,,drop = FALSE])
        } else toAdd[inWindowAdd,,drop = FALSE]
      }


    } else {    # if nr is not > 0
      added <- FALSE
      nodein <- nPairs + 1
    }
    # -----------------------------------------------------------

    # check -----------------------------------------------------
    niter <- niter + 1
    crit <- ( added == TRUE | removed == TRUE )
  } # while

  # if ( parallel ) stopCluster(cl)
  out <- list(sigma = Sigma,
              omega = Omega,
              graph = adjMat,
              loglik = critVal + penalty(adjMat, beta = beta),
              nPar = sum(adjMat)/2 + V,
              crit = critVal,
              penalty = attr(penalty, "type"),
              occam = list(TOADD = TOADD, TODROP = TODROP))
  return(out)
}
