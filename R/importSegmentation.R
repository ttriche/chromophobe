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
    if(all(levels(mcols(x)[, statecol]) %in% names(states))) {
      levels(mcols(x)[, statecol]) <- states[levels(mcols(x)[, statecol])]
    } else {
      message("The HMM has more states than you supplied, or different names:")
      stop(paste(setdiff(levels(mcols(x)[,statecol]),names(states)), coll=','))
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
  return(x)
}
