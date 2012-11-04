loadChromHMM <- function(path='.', g=NULL) { 
  if(path != '.') oldwd <- getwd()
  setwd(path)
  argv <- list()
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
            rowData=GRangesList(lapply(files, importSegmentation,loud=T,gen=g)),
            colData=DataFrame(segmentationName=names(files), bedFile=files),
            exptData=SimpleList(argv))
  if(path != '.') setwd(oldwd)
  colnames(x) <- names(files)
  return(x)
}

parseCommand <- function(command) { # {{{ parse command from webpage.html
  argv <- list()
  settings <- strsplit(strsplit(command, ':')[[1]][2], ' ')[[1]][-1]
  ## typical: "LearnModel -s 123 -r 25 -x 15000 HMM.binary HMM.model20 20 hg19"
  argv$command <- settings[1]
  if('-x' %in% settings) argv[['maxtime']] <- settings[grep('-x', settings)+1]
  if('-r' %in% settings) argv[['maxiter']] <- settings[grep('-r', settings)+1]
  if('-s' %in% settings) argv[['seed']] <- settings[grep('-s', settings)+1]
  if('-printposterior' %in% settings) argv$posterior <- TRUE
  mandatory <- settings[(length(settings)-3):length(settings)]
  names(mandatory) <- c('binaryDir','modelDir','modelStates','genome')
  for(nm in names(mandatory)) argv[[nm]] <- as.character(mandatory[nm])
  return(argv) # experimentData?
} # }}}
