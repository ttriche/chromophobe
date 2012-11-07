writeLabelStates <- function(states, file='labelStates.txt', colors=NULL) {

  ## This function is a hack so that ChromHMM 1.0.5 will relabel states nicely.
  ## Nevertheless, it's a nice hack, especially when you're making slides :-)
  labelStates <- as.data.frame(cbind(names(states), states))
  if(!is.null(colors)) message('Color support has not been added yet')
  writeDelim(labelStates, file=file, col.names=FALSE)
  message(paste('state labels are in', file, 'formatted for MakeBrowserFiles'))

}
