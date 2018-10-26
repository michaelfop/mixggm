.onAttach <- function(lib, pkg)
{
  version <- packageVersion("mixggm")
  if ( interactive() )
  { # ogre
    packageStartupMessage("
           _                           
 _ __ ___ (_)_  ____ _  __ _ _ __ ___  
| '_ ` _ \\| \\ \\/ / _` |/ _` | '_ ` _ \\    Mixtures of Gaussian
| | | | | | |>  < (_| | (_| | | | | | |   Graphical Models
|_| |_| |_|_/_/\\_\\__, |\\__, |_| |_| |_|
                 |___/ |___/              Version ", version, "\n")
  }
  else
  { packageStartupMessage("Package 'mixggm' version ", version) }

  packageStartupMessage("Type 'citation(\"mixggm\")' for citing this R package in publications.")
  invisible()
}
