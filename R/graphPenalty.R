#
#===================================== Set of graph penalty functions
#


graphPenalty <- function( penalty = c("bic", "ebic", "erdos", "power") )
  # Set the penalty to be used for graph search.
  # If the penalty is user-defined, it must be component separable
{
  if ( is.character(penalty) ) {
    penalty <- match.arg( penalty, c("bic", "ebic", "erdos", "power") )
    penFun <- switch(penalty,
                     bic = function(graph, beta = NULL)
                       # BIC type penalty
                     {
                       NN <- get("N", envir = parent.frame())
                       E <- sum(graph)/2
                       0.5*log(NN) * E
                     },
                     ebic = function(graph, beta = NULL)
                       # EBIC type penalty
                     {
                       if ( !is.null(beta) ) {
                         if ( beta > 1 | beta < 0 ) stop("The beta hyperparameter must be in the interval [0,1] for the ebic penalty function")
                       }
                       V <- ncol(graph)
                       NN <- get("N", envir = parent.frame())
                       if ( is.null(beta) ) beta <- 1
                       E <- sum(graph)/2
                       0.5*log(NN)*E + 2*beta*log(V)*E
                     },
                     erdos = function(graph, beta = NULL)
                       # Erdos-Renyi
                     {
                       if ( !is.null(beta) ) {
                         if ( beta >= 1 | beta <= 0 ) stop("The beta hyperparameter must be in the interval (0,1) for the erdos penalty function")
                       }
                       V <- ncol(graph)
                       if ( is.null(beta) ) beta <- log(V)/choose(V,2)   # expect log(V) arcs
                       TOT <- choose(V, 2)
                       E <- sum(graph)/2
                       -E*log(beta) - (TOT - E)*log(1 - beta)
                     },
                     power = function(graph, beta = NULL)
                       # power law on degree
                     {
                       NN <- get("N", envir = parent.frame())
                       if ( is.null(beta) ) beta <- 2*log(NN)
                       deg <- rowSums(graph)
                       + beta * sum( log(deg + 1) )
                     }
                     # structure = function(graph, beta = NULL)
                     #   # prior structure
                     # {
                     #   V <- ncol(graph)
                     #   if ( is.null(beta) ) {
                     #     beta <- list()
                     #     beta$coeff <- 1
                     #     beta$ref <- toeplitz( c(0, 1, rep(0, V - 2)) )
                     #   }
                     #   -beta$coeff * sum( abs(graph - beta$ref) )
                     # }
    )
    attr(penFun, "type") <- gsub("Pen", "", penalty)
  } else if ( is.function(penalty) ) {
    penFun <- penalty
    attr(penFun, "type") <- "user"
  }
  return(penFun)
}
