## realistically, if trans/emis can't be estimated, we've got a problem...
loadMethylHMM <- function(path='.',genome=NULL,states=NULL,files=NULL,para=T) { 

  require(rtracklayer)
  if(path != '.') oldwd <- getwd()
  setwd(path)
  argv <- list()
  if(is.null(files)) files <- list.files(patt='.*ethylHMM.*bed$')
  if(is.null(names(files))) names(files) <- files

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
  emis <- NULL
  trans <- NULL
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
