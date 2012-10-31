import.ChromHMM <- function(file, genome=NULL, states=NULL, statecol='name', as=c('GR','GRL')) {
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
    seqinfo(x) <- SeqinfoForBSGenome(genome)[seqlevels(x)]
  }
  names(mcols(x))[which(names(mcols(x))==statecol)] <- 'state'
  if(as == 'GRL') return(split(x, mcols(x)$state))
  else if(as == 'GR') return(x)
  else stop(paste("Don't know how to return a", as))
}

