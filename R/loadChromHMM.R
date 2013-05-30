loadChromHMM <- function(path='.',genome=NULL,states=NULL,files=NULL,para=F) { 

  require(rtracklayer)
  if(path != '.') oldwd <- getwd()
  setwd(path)
  argv <- list()
  if(is.null(files)) {
    files <- list.files(patt='.*_segments.bed$')
    names(files) <- gsub('_segments.bed$', '', files)
  }
  if(is.null(names(files))) names(files) <- files

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
    genome <- argv$genome
  }
  emis <- NULL
  trans <- NULL

  modelFile <- grep('^model.*?\\.txt$', list.files(), value=T)
  if(length(modelFile) > 0) {
    model <- importModel(modelFile, loud=TRUE)
    trans <- model$transitions
    emis <- model$emissions
    ## detect states from model, if possible and none provided
    if(is.null(states) && buildStates) states <- model$states
  } else {
    message('No model file found, omitting shared model parameters')
  } 

  if(is.null(names(files))) names(files) <- files
  if(para == TRUE) {
    require(parallel)
    segs <- mclapply(files, importSegmentation, 
                     genome=genome, states=states, loud=T)
  } else { 
    segs <- lapply(files, importSegmentation, 
                   gen=genome, states=states, loud=T)
  }
  mcol.idx <- unlist(lapply(segs, function(x) length(names(mcols(x)))))
  minMcols <- function(x, minCols=1) { # {{{
                mcols(x) <- mcols(x)[, seq_len(minCols)]
                names(mcols(x))[1] <- 'state'
                return(x)
              } # }}}
  segs <- lapply(segs, minMcols, minCols=min(mcol.idx))
  segList <- SegmentationList(GRangesList(segs), s=states)
  colDat <- DataFrame(segmentationName=names(files), bedFile=files)
  rownames(colDat) <- names(files)

  x <- new('JointSegmentation', emissions=emis, transitions=trans,
           rowData=segList, colData=colDat, states=states,
           exptData=SimpleList(argv)) 
  if(path != '.') setwd(oldwd)
  colnames(x) <- names(files)
  return(x)
}
