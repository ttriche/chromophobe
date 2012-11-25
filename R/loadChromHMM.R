loadChromHMM <- function(path='.', g=NULL, states=NULL) { 

  require(rtracklayer)
  if(path != '.') oldwd <- getwd()
  setwd(path)
  argv <- list()
  files <- list.files(patt='.*_segments.bed$')
  names(files) <- gsub('_segments.bed$', '', files)

  ## use the color-coded, labeled files if present...
  #
  # densefiles <- sapply(files, function(x) 
  #                     list.files(patt=gsub('_segments','_dense',x)))
  # if(length(densefiles) == length(files)) files <- densefiles
  #
  buildStates <- TRUE
  if(is.null(states)) {
    if('labelStates.txt' %in% list.files()) {
      states <- importLabelStates('labelStates.txt')
      buildStates <- FALSE
    }
  } else { 
    if(!is(states, 'States')) {
      states <- DataFrame(Id=names(states), Name=states)
      states <- States(states) 
      buildStates <- FALSE
    }
  }
  if(length(list.files(patt='webpage')) > 0) {
    calls <- system2('grep', c('command', list.files(patt='webpage')), stdout=T)
    argv <- parseCommand(calls)
    g <- argv$genome
  }
  modelFile <- grep('^model_?[123456789]*\\.txt$', list.files(), value=T)
  model <- importModel(modelFile, loud=TRUE)

  ## detect states from model, if possible and none provided
  if(is.null(states) && buildStates) states <- model$states

  x <- new('JointSegmentation',
            emissions=model$emissions,
            transitions=model$transitions,
            rowData=SegmentationList(lapply(files, 
                                            importSegmentation,
                                            states=states,
                                            loud=T,
                                            genome=g)),
            colData=DataFrame(segmentationName=names(files), bedFile=files),
            exptData=SimpleList(argv),
            states=states)
  if(path != '.') setwd(oldwd)
  colnames(x) <- names(files)
  return(x)
}
