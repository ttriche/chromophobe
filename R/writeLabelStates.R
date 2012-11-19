writeLabelStates <- function(states, file='labelStates.txt', colors=NULL) {

  ## This function is a hack so that ChromHMM 1.0.5 will relabel states nicely.
  ## Nevertheless, it's a nice hack, especially when you're making slides :-)
  labelStates <- as.data.frame(States)
  if(!is.null(colors) && (!'color' %in% tolower(names(labelStates)))) {
    if(length(colors) == nrow(labelStates)) labelStates$Color <- colors
  }
  writeDelim(labelStates, file=file, col.names=FALSE, row.names=TRUE)
  message(paste('state labels are in', file, 'formatted for MakeBrowserFiles'))

}
