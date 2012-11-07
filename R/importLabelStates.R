importLabelStates <- function(file) {
  require(utils)
  s <- read.delim(file, header=F)
  states <- s[,2]
  names(states) <- s[,1]
  return(states)
}
