mixggmStartupMessage <- function()
{
  msg <- c(paste0(
"           _                           
 _ __ ___ (_)_  ____ _  __ _ _ __ ___  
| '_ ` _ \\| \\ \\/ / _` |/ _` | '_ ` _ \\    Mixtures of Gaussian
| | | | | | |>  < (_| | (_| | | | | | |   Graphical Models
|_| |_| |_|_/_/\\_\\__, |\\__, |_| |_| |_|
                 |___/ |___/              Version ",
  packageVersion("mixggm")),
  "\nType 'citation(\"mixggm\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- mixggmStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'mixggm' version", packageVersion("mixggm"))
  packageStartupMessage(msg)      
  invisible()
}