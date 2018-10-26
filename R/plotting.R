#
#======================== Plotting functionalities
#
#

plot.mixGGM <- function(x, what = c("graph", "classification", "adjacency", "common"),
                        layout = c("circle", "random"), colors = NULL, symb = NULL, dimens = NULL, ...)
{
  what <- match.arg( what, c("graph", "classification", "adjacency", "common") )
  layout <- match.arg( layout, c("circle", "random") )

  if ( is.null(colors) ) {
    colors <- c("#AA4488", "#4477AA", "#AA7744", "#AAAA44", "#AA4455", "#44AA77", "#44AAAA",
                "#771155", "#114477", "#774411", "#777711", "#771122", "#117744", "#117777",
                "#CC99BB", "#77AADD", "#DDAA77", "#DDDD77", "#DD7788", "#88CCAA", "#77CCCC" )
  }
  if ( is.null(symb) ) {
    symb <- c(16, 15, 17, 1, 3, 0, 8, 2, 4,  7,  5,  9,  6, 10, 11, 18, 12, 13, 14, 19, 20, 21)
  }
  graph <- x$graph
  sigma <- x$parameters$sigma
  omega <- x$parameters$omega
  mu <- x$parameters$mu
  K <- x$K
  varnames <- if ( !is.null(colnames(sigma[,,1])) ) colnames(sigma[,,1]) else paste0("V", 1:ncol(sigma[,,1]))
  V <- length(varnames)
  op <- par()

  # graphs ----------------------------------------------------------
  if ( what == "graph" ) {
    #
    if ( K > 3 ) {
      par( mfrow = c(2,3), mar = rep(0.2,4) )
    } else if ( K > 1 ) {
      par( mfrow = c(ifelse(K == 2, 1, 2),2), mar = rep(0.2,4) )
    }
    tmp <- network::network(graph[,,1])
    coord <- if ( layout == "circle" ) network::network.layout.circle(tmp) else NULL
    r <- if ( x$model == "covariance" ) cov2cor( sigma[,,1] ) else cov2cor( omega[,,1] )
    network::plot.network( tmp, coord = coord, label = varnames, jitter = FALSE, edge.lwd = r*10,
                           vertex.col = adjustcolor(colors[1], 0.8), edge.col = adjustcolor("gray60", 0.6),
                           arrowhead.cex = 1.5, vertex.cex = 1.5,
                           label.pos = 3, label.cex = 1, usearrows = x$model == "covariance", ... )
    if ( K > 1 ) {
      for ( h in 2:K ) {
        tmp <- network::network(graph[,,h])
        r <- if ( x$model == "covariance" ) cov2cor( sigma[,,h] ) else cov2cor( omega[,,h] )
        network::plot.network( tmp, coord = coord, label = varnames, jitter = FALSE, edge.lwd = r*10,
                               vertex.col = adjustcolor(colors[h], 0.8), edge.col = adjustcolor("gray60", 0.6),
                               arrowhead.cex = 1.5, vertex.cex = 1.5, label.pos = 3, label.cex = 1,
                               usearrows = x$model == "covariance", ...)
      }
    }
    on.exit( par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1) )

    # adjacency -------------------------------------------------------
  } else if ( what == "adjacency" ) {
    if ( K > 4 ) {
      par( mfrow = c(2,4), mar = rep(0.9,4) )
    } else if ( K > 1 ) {
      par( mfrow = c(ifelse(K == 2, 1, 2),2), mar = rep(0.9,4) )
    }

    for ( h in 1:K ) {
      image( 1:V, 1:V, graph[,V:1,h], col = c("white", colors[h]),
             axes = FALSE, xlab = "", ylab = "" )
      abline(V+1,-1, lwd = 2.5, col = "gray60")
      abline(v = 1:(V-1) + 0.5, col = "gray60", lwd = 0.7)
      abline(h = 1:(V-1) + 0.5, col = "gray60", lwd = 0.7)
      abline(v = c(0,V) + 0.5, col = "gray60")
      abline(h = c(0,V) + 0.5, col = "gray60")
    }
    on.exit( par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1) )

    # classification -------------------------------------------------------
  } else if ( what == "classification" ){
    #
    class(x) <- "Mclust"
    x$parameters$variance$sigma <- sigma
    x$parameters$mean <- mu
    x$parameters$variance$cholsigma <-
      array( apply(sigma, 3, chol), c(V,V,K) )
    mclust::plot.Mclust(x, what = "classification", dimens = dimens,
                        symbols = symb, colors = adjustcolor(colors, 0.9))

    # common edges --------------------------------------------------
  } else if ( what == "common" ) {
    if ( K < 2 ) stop("Need at least 2 clusters")
    g <- graph[,,1]
    for ( h in 2:K )  g <- g * graph[,,h]
    diag(g) <- 0
    tmp <- network::network(g)
    coord <- network::network.layout.circle(tmp)
    network::plot.network( tmp, coord = coord, label = varnames, jitter = FALSE,
                           vertex.col = "gray10", edge.col = adjustcolor("gray60", 0.6),
                           vertex.cex = 1.5, label.pos = 3, label.cex = 1, usearrows = x$model == "covariance")
  }
}



plot.fitGGM <- function(x, what = c("graph", "adjacency"),
                        layout = c("circle", "random"), ...)
{

  what <- match.arg( what, c("graph", "adjacency") )
  layout <- match.arg( layout, c("circle", "random") )
  class(x) <- "mixGGM"
  varnames <- dimnames(x$sigma)
  x$K <- 1
  x$parameters$mu <- 0
  x$parameters$sigma <- array(x$sigma, dim = c(x$V, x$V,1))
  x$parameters$omega <- array(x$omega, dim = c(x$V, x$V,1))
  x$graph <- array(x$graph, dim = c(x$V,x$V,1))
  dimnames(x$parameters$sigma) <- dimnames(x$graph) <- varnames
  plot.mixGGM(x, what = what, layout = layout, colors = "gray20")
}
