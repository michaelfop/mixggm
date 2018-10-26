#
#======================== EM algorithm for mixture of Gaussian graphical models
#
#

emMixGGM <- function( data, K, search = c("step-forw", "step-back","ga"),
                      model = c("covariance", "concentration"),
                      penalty = graphPenalty(),
                      beta = NULL,
                      regularize = FALSE,
                      regHyperPar = NULL,
                      ctrlEm = ctrlEM(),
                      ctrlStep = ctrlSTEP(),
                      ctrlGa = ctrlGA(),
                      ctrlIcf = ctrlICF(),
                      hcInit = NULL,
                      parallel = FALSE,
                      verbose = FALSE )
{
  varnames <- colnames(data)
  N <- nrow(data)
  V <- ncol(data)
  search <- match.arg( search,  c("step-forw", "step-back", "ga") )
  data <- data.matrix(data)

  # regularization? -------------------------------------------------
  if ( regularize ) {
    scaleType <- attr(regularize, "scaleType")
    scaleType <- if ( is.null(scaleType) ) "full" else scaleType
    if ( is.null(regHyperPar) ) regHyperPar <- ctrlREG(data, K, scaleType = scaleType)
    class(regHyperPar) <- "EM"
  }
  #------------------------------------------------------------------

  ##### TODO: NEED TO BE CHECKED
  # start parallel computations -------------------------------------
  if ( !inherits(parallel, "cluster") ) {    # if parallel not already started then
    if ( parallel | is.numeric(parallel) ) {
      parallel <- GA::startParallel(parallel)
      class(parallel) <- "cluster"
    }
    on.exit( if ( parallel ) parallel::stopCluster(attr(parallel, "cluster")) )
  } # else cluster will be closed in 'mixGGraph' function
  #------------------------------------------------------------------

  # EM initialization -----------------------------------------------
  if ( is.null(hcInit) ) {
    # function not called from mixGGraph
    if ( !is.null( ctrlEm$subset ) ) {            # initialize from subset of data
      hcInit <- mclust::hc(data[ctrlEm$subset,], modelName = "VVV", use = "VARS")
      z <- mclust::unmap( mclust::hclass(hcInit, K) )
      ms <- mclust::mstep(modelName = "VVV", z = z, data = data[ctrlEm$subset,])
      es <- do.call( "estep", c(list(data = data), ms) )
      z <- es$z
    } else {
      hcInit <- mclust::hc(data, modelName = "VVV", use = "VARS")
      z <- mclust::unmap( mclust::hclass(hcInit, K) )
    }
  } else {
    # function called from mixGGraph
    if ( !is.null(ctrlEm$subset) ) {             # initialize from subset of data
      z <- mclust::unmap( mclust::hclass(hcInit, K) )
      ms <- mclust::mstep(modelName = "VVV", z = z, data = data[ctrlEm$subset,])
      es <- do.call( "estep", c(list(data = data), ms) )
      z <- es$z
    } else {
      z <- mclust::unmap( mclust::hclass(hcInit, K) )
    }
  }
  pro <- colMeans(z)

  # automatic regularization if Nk < M
  if ( any( table(mclust::map(z)) <= V ) ) {
    if ( !regularize ) regularize <- TRUE
    if ( is.null(regHyperPar) ) {
      regHyperPar <- ctrlREG(data, K)
      class(regHyperPar) <- "EM"
    }
  }
  #-------------------------------------------------------------------------

  # EM and starting parameters ---------------------------------------------
  tol <- ctrlEm$tol
  itMax <- ctrlEm$maxiter
  #
  iter <- 0
  loglikPrev <- -.Machine$integer.max/2
  loglikPenPrev <- -.Machine$integer.max/2
  err <- .Machine$double.xmax/2
  if ( is.null(ctrlStep$start) ) {
    firstIt <- TRUE
    ctrlStep$start <- 0
  } else firstIt <- FALSE    # if start matrix is given
  Start <- array( ctrlStep$start, c(V,V,K) )
  occam <- rep( list(list(TOADD = NULL, TODROP = NULL)), K )
  critOut <- rep(NA, K)
  crit <- TRUE
  #-------------------------------------------------------------------------

  # EM ---------------------------------------------------------------------
  if ( regularize ) {
    Scale <- regHyperPar$scale
    cholScale <- chol(Scale)
  }
  sigma <- omega <- graph <- array( NA, c(V,V,K) )
  dimnames(sigma) <- dimnames(omega) <- dimnames(graph) <- list(varnames, varnames)

  while ( crit ) {

    #### M step ............................................
    temp <- mclust::covw(data, z, normalize = FALSE)
    dimnames(temp$S) <- list(varnames, varnames)
    mu <- temp$mean
    pro <- colMeans(z)

    for ( j in 1:K ) {
      if ( firstIt ) START <- NULL else {
        START <- structure(Start[,,j], critOut = critOut[j], sigma = sigma[,,j])
      }
      Nj <- pro[j]*N

      if ( regularize ) {
        Sj <- Scale + Nj*temp$S[,,j]
        Sj <- Sj / ( Nj + V + 1 + regHyperPar$psi )
      } else {
        Sj <- temp$S[,,j]
      }

      out <- switch(search,
                    "step-forw" = try(
                      searchGGMStepwise_f(S = Sj, N = N, model = model, pro = Nj/N, start = START,
                                                 penalty = penalty, beta = beta,
                                                 regularize = regularize, regHyperPar = regHyperPar,
                                                 ctrlStep = ctrlStep, ctrlIcf = ctrlIcf, parallel = parallel,
                                                 verbose = FALSE,
                                                 occam = occam[[j]]),
                      silent = ctrlEm$printMsg),
                    #
                    "step-back" = try(
                      searchGGMStepwise_b(S = Sj, N = N, model = model, pro = Nj/N, start = START,
                                                 penalty = penalty, beta = beta,
                                                 regularize = regularize, regHyperPar = regHyperPar,
                                                 ctrlStep = ctrlStep, ctrlIcf = ctrlIcf, parallel = parallel,
                                                 verbose = FALSE,
                                                 occam = occam[[j]]),
                      silent = !ctrlEm$printMsg),
                    #
                    "ga" = try(
                      searchGGMGA(S = Sj, N = N, model = model, pro = Nj/N, start = START,
                                         penalty = penalty, beta = beta,
                                         regularize = regularize, regHyperPar = regHyperPar,
                                         ctrlGa = ctrlGa, ctrlIcf = ctrlIcf, parallel = parallel,
                                         verbose = FALSE),
                      silent = !ctrlEm$printMsg)
      )

      if ( class(out) == "try-error" ) {
        exitLoop <- TRUE
        break
      } else exitLoop <- FALSE

      sigma[,,j] <- out$sigma
      omega[,,j] <- out$omega
      graph[,,j] <- out$graph
      occam[[j]] <- out$occam   # keep the edges in the window for the next EM iteration
      critOut[j] <- out$crit

    } # for j in 1:K
    if ( exitLoop ) break
    ###.....................................................

    #### E step ............................................
    e <- estepmggm(data, t(mu), sigma, pro)
    z <- e$z
    loglik <- e$loglik
    #
    ### OLD VERSION
    # cholSigma <- array( apply(sigma, 3, chol), c(V,V,K) )
    # e <- estepVVV(x, parameters = list(mean = mu, pro = pro,
    #                                    variance = list(cholsigma = cholSigma))
    #               )
    ###.....................................................

    ### loglikelihood.......................................
    llk <- loglik             # BIC is computed using MAP estimates in loglikelihood
    #
    if ( regularize ) {
      nu <- regHyperPar$psi
      logPrior <- rep(NA, K)
      for ( j in 1:K ) {
        # inverse Wishart log density  --- MCMCpack
        multGamma <- sum( lgamma((nu + 1 - 1:V)/2) )
        lDenom <- multGamma + 0.5*nu*V*log(2) + 0.25*V*(V - 1)*log(pi)
        cholS <- chol(sigma[,,j])
        cholW <- cholScale
        halflogdetS <- sum( log(diag(cholS)) )
        halflogdetW <- sum( log(diag(cholW)) )
        invS <- chol2inv( cholS )
        exptrace <- sum(Scale * invS)
        lNum <- nu*halflogdetW - (nu + V + 1) * halflogdetS - 0.5 * exptrace
        logPrior[j] <- lNum - lDenom
      }
      loglik <- loglik + sum(logPrior)
    }

    nCov <- sum(graph)/2      # number of covariance parameters
    penVal <- sum( vapply( 1:K, function(k) penalty(graph[,,k], beta = beta), numeric(1) ) )
    loglikPen <- loglik - penVal
    ###.....................................................

    # update the Start......................................
    # temp <- ifelse(sigma != 0, 1, 0)
    # temp <- apply( temp, 3, function(a) replace(a, cbind(1:M, 1:M), 0) )
    Start <- graph
    firstIt <- FALSE
    #.......................................................

    #### check
    err <- abs(loglikPen - loglikPenPrev) / (1 + abs(loglikPen))
    loglikPenPrev <- loglikPen
    iter <- iter + 1

    # use loading bar
    crit <- ( err > tol & iter < itMax )
  } # while
  #---------------------------------------------------------------------------

  ### extra M step ...........................................
  # we are implementing a ME algorithm so we need an extra M step!
  if ( exitLoop ) {
    pro[] <- mu[] <- sigma[] <- omega[] <- graph[] <- z[] <- NA
    loglik <- llk <- loglikPen <- nCov <- totPar <- NA
  } else {
    temp <- mclust::covw(data, z, normalize = FALSE)
    dimnames(temp$S) <- list(varnames, varnames)
    mu <- temp$mean
    #
    for ( j in 1:K ) {
      Nj <- sum(z[,j])
      if ( regularize ) {
        Sj <- Scale + Nj*temp$S[,,j]
        Sj <- Sj / ( Nj + V + 1 + regHyperPar$psi )
      } else {
        Sj <- temp$S[,,j]
      }
      tmp <- fitGGM( data = NULL, graph = graph[,,j], S = Sj, model = model, N = Nj,
                     ctrlIcf = ctrlIcf, regularize = regularize, regHyperPar = regHyperPar )
      sigma[,,j] <- tmp$sigma
      omega[,,j] <- tmp$omega
    }
    #
    pro <- colMeans(z)
  }
  ###.........................................................

  dimnames(mu) <- list(varnames, 1:K)
  dimnames(sigma) <- dimnames(omega) <- dimnames(graph) <- list(varnames, varnames)
  totPar <- V*K + V*K + nCov + (K-1)
  out <- list( parameters = list(tau = pro, mu = mu, sigma = sigma, omega = omega), graph = graph,
               N = N, V = V, K = K, loglik = loglik, loglikPen = loglikPen, loglikReg = llk,
               nPar = c(depPar = sum(nCov), totPar = totPar), z = z, classification = mclust::map(z),
               penalty = attr(penalty, "type") )
  attr(out, "info") <- list(iter = iter, control = ctrlEm[1:2], search = search)
  return(out)
}
