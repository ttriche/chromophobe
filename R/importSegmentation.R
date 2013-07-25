## take a segmentation (typically _dense.bed or _segments.bed) and process it 
importSegmentation <- function(file, 
                               genome=NULL,
                               states=NULL, 
                               statecol='name', 
                               loud=FALSE,
                               addGaps=FALSE,
                               nonGapState=NULL,
                               gapState='OTHER') {
  require(rtracklayer)
  if(loud) message(paste('Importing segmentation from', file, '...'))
  x <- import.bed(file, asRangedData=FALSE)

  ## for MethylHMMs and DNAseI peaks/footprints
  if(addGaps==TRUE) {
    if(!is.null(nonGapState)) {
      mcols(x)[, statecol] <- nonGapState
      mcols(x)[, 'type'] <- nonGapState
    } else {
      mcols(x)[, 'type'] <- mcols(x)[, statecol]
    }
    gapsX <- gaps(x)
    gapsX <- gapsX[ which(strand(gapsX) == '*') ]
    mcols(gapsX)[, statecol] <- gapState
    mcols(gapsX)[, 'type'] <- gapState
    mcols(gapsX)[, 'score'] <- 0
    keep.cols <- c('score', statecol, 'type') 
    mcols(x)[,statecol] <- mcols(x)$type
    mcols(x) <- mcols(x)[, keep.cols]
    mcols(gapsX) <- mcols(gapsX)[, keep.cols]
    merged <- sort(c(x, gapsX))
    mcols(merged)$type <- as.factor(mcols(merged)$type)
    mcols(merged)$score <- as.numeric(mcols(merged)$score)
    x <- merged
  }

  mcols(x)[, statecol] <- as.factor(mcols(x)[, statecol])
  if(!is.null(states)) {
    stopifnot(class(states) == 'States')
    if(all(levels(mcols(x)[, statecol]) %in% stateIds(states))) {
      lvls <- stateNames(states)[levels(mcols(x)[,statecol])]
      levels(mcols(x)[,statecol]) <- lvls
    } else {
      message("The HMM has more states than you supplied, or different names:")
      message('Supplied states: ', paste(stateIds(states), collapse=', '))
      message('Supplied state names: ',paste(stateNames(states), collapse=', '))
      message('Model states: ',paste(levels(mcols(x)[,statecol]),collapse=', '))
      stop()
    }
  }
  if(!is.null(genome)) {
    genome(x) <- genome
    seqinfo(x) <- rtracklayer::SeqinfoForBSGenome(genome)[seqlevels(x)]
  } else {
    message('You did not specify a genome -- this could cause trouble later on')
  }
  names(mcols(x))[which(names(mcols(x))==statecol)] <- 'state'
  return(Segmentation(x, states))
}
