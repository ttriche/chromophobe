import.ChromHMM <- function(file, genome=NULL, states=NULL, as=c('GR','GRL')) {
  x <- import.bed(file, asRangedData=FALSE)
  mcols(x)$name <- as.factor(mcols(x)$name)
  if(!is.null(states) && all(levels(mcols(x)$name) %in% names(states))) {
    levels(mcols(HMM)$name) <- states[levels(mcols(HMM)$name)]
  }
  if(!is.null(genome)) {
    genome(x) <- genome
    seqinfo(x) <- SeqinfoForBSGenome(genome)[seqlevels(x)]
  }
  names(mcols(x))[ which(names(mcols(x)) == 'name') ] <- 'state'
  if(as == 'GRL') return(split(x, mcols(x)$state))
  else if(as == 'GR') return(x)
  else stop(paste("Don't know how to return a", as))
}

