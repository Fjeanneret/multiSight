multiSightStartupMessage <- function()
{
  msg <- "Please, find comprehensive tutorial of multiSight functions within vignette by browseVignettes('multiSight')\nThe graphical interface is described in README.md, see https://github.com/Fjeanneret/multiSight/\nFor any questions or issues, you can find help in https://github.com/Fjeanneret/multiSight/issues."
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- multiSightStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'multiSight' version", packageVersion("multiSight"))
  packageStartupMessage(msg)
  invisible()
}
