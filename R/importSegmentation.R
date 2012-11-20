## take a segmentation (typically _dense.bed or _segments.bed) and process it 
importSegmentation <- function(file, 
                               genome=NULL,
                               states=NULL, 
                               statecol='name', 
                               loud=FALSE) {
  if(loud) message(paste('Importing segmentation from', file, '...'))
  x <- import.bed(file, asRangedData=FALSE)
  mcols(x)[, statecol] <- as.factor(mcols(x)[, statecol])
  if(!is.null(states)) {
    stopifnot(class(states) == 'States')
    if(all(levels(mcols(x)[, statecol]) %in% stateIds(states))) {
      lvls <- stateNames(states)[levels(mcols(x)[,statecol])]
      levels(mcols(x)[,statecol]) <- lvls
    } else {
      message("The HMM has more states than you supplied, or different names:")
      browser()
      message('Supplied states: ', paste(stateIds(states), collapse=', '))
      message('Supplied state names: ',paste(stateNames(states), collapse=', '))
      message('Model states: ',paste(levels(mcols(x)[,statecol]),collapse=', '))
      stop()
    }
  }
  if(!is.null(genome)) {
    genome(x) <- genome
    require(rtracklayer)
    seqinfo(x) <- rtracklayer::SeqinfoForBSGenome(genome)[seqlevels(x)]
  } else {
    message('You did not specify a genome -- this could cause trouble later on')
  }
  names(mcols(x))[which(names(mcols(x))==statecol)] <- 'state'
  return(as(x, 'Segmentation'))
}
