importLabelStates <- function(file) {
  require(utils)
  s <- read.delim(file, header=F)
  states <- s[,2]
  names(states) <- s[,1]
  if(ncol(s) > 2) {
    colrs <- s[,3]
  } else {
    colrs <- rainbow(nrow(s))
  }
  DF <- DataFrame(Id=s[,1], Name=s[,2], Color=colrs)
  rownames(DF) <- s[,1]
  return(States(DF))
}
