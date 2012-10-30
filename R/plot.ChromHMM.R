## plot genome or site % occupancy for a given HMM's states
##
## FIXME: handle WIG tracks and other things with scores/depths
##
plot.ChromHMM <- function(HMMs, GR=NULL, asReads=FALSE, ...) {

  require(ggplot2)
  if(is.null(GR)) what <- 'genomic base'
  else what <- paste(as.character(match.call()["GR"]), 'site')
  title <- paste(what, 'occupancy by HMM state')
  if(is.list(HMMs)) {
    if(is.list(GR) || is(GR, 'GRangesList')) {
      message('Lists of comparisons are not yet supported, but they should be!')
    }
    byState <- data.frame(do.call(cbind, 
                                  lapply(HMMs, 
                                         basesByState, 
                                         GR=GR, 
                                         asReads=asReads)))  
  } else { 
    byState <- data.frame(HMM=basesByState(HMMs, GR=GR, asReads=asReads))
  }
  byState$state <- factor(rownames(byState))

  require(reshape2)
  byState <- melt(byState, id.vars='state')
  names(byState) <- gsub('^variable$', 'cell', names(byState))
  names(byState) <- gsub('^value$', 'bases', names(byState))
  ggplot(byState, aes(y=bases, x=factor(cell), fill=state)) + 
    geom_bar(position="fill") + 
    theme(plot.title=element_text(face="bold", size=14),
          axis.title.y=element_text(face="bold", colour="#990000", size=14),
          axis.text.y=element_text(size=12), 
          axis.title.x=element_text(face="bold", colour="#990000", size=14),
          axis.text.x=element_text(angle=45, hjust=1, size=12)) + 
    ylab(paste0('Cumulative fraction of ', what, 's occupied')) +
    xlab('Cell type') + 
    ggtitle(title)
}

basesByState <- function(x, GR=NULL, asReads=F, wig=F, statecol='state', ...) {
  stopifnot(is(x, 'GenomicRanges'))
  if(wig == TRUE) stop('Wig files and pileups are not yet supported...')
  if(!is.null(GR)) {
    if(asReads==TRUE) { ## don't disjoin and count bases if handed reads
      byState <- split(x, mcols(x)[,statecol])
      sort(unlist(lapply(byState, function(xx) sum(countOverlaps(xx, GR,...)))))
    } else { 
      ol <- findOverlaps(x, GR)
      subx <- x[queryHits(ol)]
      subGR <- GR[subjectHits(ol)]
      disjoint <- subsetByOverlaps(disjoin(c(subx,subGR,ignore.mcols=T)), subGR)
      mcols(disjoint)[,statecol] <- rep('', length(disjoint))
      byState <- split(subx, mcols(subx)[,statecol])
      for(state in names(byState)) {
        overlapping <- queryHits(findOverlaps(disjoint, byState[[state]]))
        if(length(overlapping)>0) mcols(disjoint[overlapping])[,statecol]<-state
      }
      mcols(disjoint)[,statecol] <- as.factor(mcols(disjoint)[,statecol])
      byState <- split(disjoint, mcols(disjoint)[,statecol])
      sort(unlist(lapply(byState, function(xx) sum(as.numeric(width(xx))))))
    }
  } else { 
    byState <- split(x, mcols(x)[,statecol])
    sort(unlist(lapply(byState, function(xx) sum(as.numeric(width(xx))))))
  }
}
