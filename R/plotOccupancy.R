plotOccupancy <- function(x, dropQuiescent=FALSE, stateColours=NULL) { 
  require(grid)
  require(ggplot2)
  title <- 'Chromatin state occupancy'
  x[,'state'] <- factor(rownames(x))
  if(dropQuiescent==TRUE) {
    x <- x[grep('quiescent', invert=TRUE, tolower(rownames(x))),]
    title <- paste(title, '(exclusive of quiescent states)')
  }
  colorscale <- NULL
  if(!is.null(stateColours)) colorscale <- stateColours
  else if(!is.null(stateColors(x))) colorscale <- stateColors(x)
  require(reshape2)
  byState <- melt(as.data.frame(x), id.vars='state')
  names(byState) <- gsub('^variable$', 'cell', names(byState))
  names(byState) <- gsub('^value$', 'fraction', names(byState))
  p <- ggplot(byState, aes(y=fraction, x=factor(cell), fill=state)) + 
         geom_bar(position="fill") + 
         theme(plot.title=element_text(face="bold",size=14),
               panel.background=element_blank(), 
               panel.margin = unit(0, "lines"), 
               axis.title.x=element_text(face="bold",color="#990000",size=14),
               axis.title.y=element_text(face="bold",color="#990000",size=14),
               axis.text.x=element_text(angle=45, hjust=1, vjust=1,
                                        size=10, color="#990000"), 
               axis.text.y=element_blank(),
               axis.ticks.x=element_blank(),
               axis.ticks.y=element_blank()) +
         ylab(paste0('Fraction of sites occupied')) + 
         xlab('') + 
         ggtitle(title)
  if(!is.null(colorscale)) p <- p + scale_fill_manual(values=colorscale)
  return(p)  
}

