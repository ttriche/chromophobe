writeColorStates <- function(obj, file='colorStates.txt', colors=NULL) {

  colorStates <- as.data.frame(states(obj))
  stopifnot('color' %in% tolower(names(colorStates)))
  colorStates <- colorStates[, c('Id','Color')]
  colorStates$Id <- gsub('^E', '', colorStates$Id)
  colorStates$Color <- apply(col2rgb(colorStates$Color), 2, paste, collapse=',')
  writeDelim(colorStates, file=file, col.names=FALSE, row.names=FALSE, sep="\t")
  message(paste('state colors are in', file, 'formatted for MakeBrowserFiles'))

}
