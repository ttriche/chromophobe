import.ChromHMM <- function(file, genome=NULL, states=NULL, as=c('GR','GRL')) {
  x <- import.bed(file, asRangedData=FALSE)
  mcols(x)$name <- as.factor(mcols(x)$name)
  if(!is.null(states)) {
    if(all(levels(mcols(x)$name) %in% names(states))) {
      levels(mcols(x)$name) <- states[levels(mcols(x)$name)]
    } else {
      stop("The HMM has more states than your list, or different state names.")
    }
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

