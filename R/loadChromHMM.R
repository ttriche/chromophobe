loadChromHMM <- function(path='.', g=NULL, states=NULL) { 

  if(path != '.') oldwd <- getwd()
  setwd(path)
  argv <- list()
  if(grepl('labelStates.txt', list.files()) && is.null(states)) {
    states <- importLabelStates('labelStates.txt')
  }
  if(length(list.files(patt='webpage')) > 0) {
    calls <- system2('grep', c('command', list.files(patt='webpage')), stdout=T)
    argv <- parseCommand(calls)
    g <- argv$genome
  }
  modelFile <- grep('^model_.*\\.txt$', list.files(), value=T)
  model <- importModel(modelFile, loud=TRUE)
  files <- list.files(patt='.*_segments.bed$')
  names(files) <- gsub('_segments.bed$', '', files)
  x <- new('JointSegmentation',
            emissions=model$emissions,
            transitions=model$transitions,
            rowData=GRangesList(lapply(files, 
                                       importSegmentation,
                                       states=states,
                                       loud=T,
                                       gen=g)),
            colData=DataFrame(segmentationName=names(files), bedFile=files),
            exptData=SimpleList(argv))
  if(path != '.') setwd(oldwd)
  colnames(x) <- names(files)
  return(x)

}
