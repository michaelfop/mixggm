#
#======================== Plotting functionalities
#
#

plot.mixGGM <- function(x, 
                        what = c("graphs", "classification", "adjacency", "common"),
                        layout = c("circle", "random"), 
                        col = mclust::mclust.options("classPlotColors"), 
                        pch = mclust::mclust.options("classPlotSymbols"), 
                        dimens = NULL, 
                        ...)
{
  what   <- match.arg(what, choices = eval(formals(plot.mixGGM)$what),
                            several.ok = TRUE)
  layout <- match.arg(layout, choices = eval(formals(plot.mixGGM)$layout))

  args  <- list(...)
  graph <- x$graph
  sigma <- x$parameters$sigma
  omega <- x$parameters$omega
  mu <- x$parameters$mu
  K <- x$K
  varnames <- if( !is.null(colnames(sigma[,,1])) ) 
                colnames(sigma[,,1]) else paste0("V", 1:ncol(sigma[,,1]))
  V <- length(varnames)
  
  # graph ----------------------------------------------------------
  plot_graphs <- function(...)
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if( K > 3 ) 
    { par( mfrow = c(2,3), mar = rep(1,4) )
    } else 
    if ( K > 1 ) 
    { par( mfrow = c(ifelse(K == 2, 1, 2),2), mar = rep(1,4) )
    }
    tmp <- network::network(graph[,,1])
    if(is.null(args$coord))
    { coord <- if(layout == "circle") 
                  network::network.layout.circle(tmp) else NULL 
    } else 
    { coord <- args$coord
      args$coord <- NULL 
    }
    r <- if(x$model == "covariance") cov2cor(sigma[,,1]) else cov2cor(omega[,,1])
    do.call("plot.network",
            c(list(tmp, label = varnames, 
                   coord = coord, # mode = "circle",
                   jitter = FALSE, 
                   vertex.col = adjustcolor(col[1], 0.8), 
                   vertex.border = col[1],
                   edge.col = adjustcolor("gray60", 0.6),
                   edge.lwd = r*10,
                   arrowhead.cex = 1.5, vertex.cex = 1.5,
                   label.pos = 3, label.cex = 1, 
                   usearrows = x$model == "covariance"),
              args))
    for(k in seq_len(K)[-1])
    {
      tmp <- network::network(graph[,,k])
      r <- if (x$model == "covariance") cov2cor(sigma[,,k]) else cov2cor(omega[,,k])
      do.call("plot.network",
              c(list(tmp, label = varnames, 
                     coord = coord, jitter = FALSE, 
                     vertex.col = adjustcolor(col[k], 0.8), 
                     vertex.border = col[k],
                     edge.col = adjustcolor("gray60", 0.6),
                     edge.lwd = r*10,
                     arrowhead.cex = 1.5, vertex.cex = 1.5, 
                     label.pos = 3, label.cex = 1,
                     usearrows = x$model == "covariance"),
                args))
    }
  }
 
  # adjacency -------------------------------------------------------
  plot_adjacency <- function(...)
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if ( K > 4 ) 
    { par( mfrow = c(2,4), mar = rep(1,4), pty = "s")
    } else 
    if ( K > 1 ) 
    {
      par( mfrow = c(ifelse(K == 2, 1, 2),2), mar = rep(1,4), pty = "s")
    }

    for(k in 1:K) 
    {
      image(1:V, 1:V, graph[,V:1,k], col = c("white", col[k]),
            axes = FALSE, xlab = "", ylab = "" )
      abline(V+1, -1, lwd = 2.5, col = "gray60")
      abline(v = 1:(V-1) + 0.5, col = "gray60", lwd = 0.7)
      abline(h = 1:(V-1) + 0.5, col = "gray60", lwd = 0.7)
      abline(v = c(0,V) + 0.5, col = "gray60")
      abline(h = c(0,V) + 0.5, col = "gray60")
    }
  }

  # classification -------------------------------------------------------
  plot_classification <- function(...)
  {
    class(x) <- "Mclust"
    x$parameters$variance$sigma <- sigma
    x$parameters$mean <- mu
    x$parameters$variance$cholsigma <-
      array( apply(sigma, 3, chol), c(V,V,K) )
    mclust::plot.Mclust(x, what = "classification", dimens = dimens,
                        symbols = pch, colors = adjustcolor(col, 0.9))
  }
  # common edges --------------------------------------------------
  plot_common <- function(...)
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1,1), mar = rep(1,4))
    if ( K < 2 ) stop("Need at least 2 clusters")
    g <- graph[,,1]
    for (k in seq_len(K)[-1])  g <- g * graph[,,k]
    diag(g) <- 0
    tmp <- network::network(g)
    coord <- network::network.layout.circle(tmp)
    network::plot.network(tmp, coord = coord, label = varnames, 
                          jitter = FALSE,
                          vertex.col = "gray10", 
                          edge.col = adjustcolor("gray60", 0.6),
                          vertex.cex = 1.5, 
                          label.pos = 3, label.cex = 1, 
                          usearrows = x$model == "covariance")
  }
  
  if(interactive() & length(what) > 1)
  { 
    title <- "Mixture of GGMs plots:"
    # present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
    while(choice != 0)
    { if(what[choice] == "graphs")         plot_graphs(...)
      if(what[choice] == "classification") plot_classification(...)
      if(what[choice] == "adjacency")      plot_adjacency(...)
      if(what[choice] == "common")         plot_common(...)
      # re-present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
    }
  } else 
  { 
    if(any(what == "graphs"))         plot_graphs(...)
    if(any(what == "classification")) plot_classification(...) 
    if(any(what == "adjacency"))      plot_adjacency(...) 
    if(any(what == "common"))         plot_common(...) 
  }
  invisible()
}



plot.fitGGM <- function(x, what = c("graphs", "adjacency"),
                        layout = c("circle", "random"), ...)
{

  what <- match.arg( what, c("graphs", "adjacency") )
  layout <- match.arg( layout, c("circle", "random") )
  class(x) <- "mixGGM"
  varnames <- dimnames(x$sigma)
  x$K <- 1
  x$parameters$mu <- 0
  x$parameters$sigma <- array(x$sigma, dim = c(x$V, x$V,1))
  x$parameters$omega <- array(x$omega, dim = c(x$V, x$V,1))
  x$graph <- array(x$graph, dim = c(x$V,x$V,1))
  dimnames(x$parameters$sigma) <- dimnames(x$graph) <- varnames
  plot.mixGGM(x, what = what, layout = layout, col = "gray20")
}
