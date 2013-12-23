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

  ## have to explicitly namespace this to avoid ggbio collision
  p <- ggplot2::ggplot(byState, 
                       aes(factor(cell), fraction, fill=state, geom="bar")) +
         ggplot2::geom_bar(position="fill") + 
         ggplot2::ylab('Fraction of sites occupied') + 
         ggplot2::xlab('') + 
         ggplot2::ggtitle(title)

  ## theme separately
  p <- p + ggplot2::theme(plot.title=element_text(face="bold",size=14),
                          panel.background=element_blank(), 
                          panel.margin = unit(0, "lines"), 
                          axis.title.x=element_text(face="bold",
                                                    color="#990000",
                                                    size=14),
                          axis.title.y=element_text(face="bold",
                                                    color="#990000",
                                                    size=14),
                          axis.text.x=element_text(angle=45, 
                                                   hjust=1, 
                                                   vjust=1,
                                                   size=10, 
                                                   color="#990000"), 
                          axis.text.y=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.ticks.y=element_blank())

  ## add state colors
  if(!is.null(colorscale)) {
    p <- p + ggplot2::scale_fill_manual(values=colorscale)
  }

  return(p)  
}

## getOccBySign -- in parallel
getOccBySign <- function(DMRs, y) { # {{{
  require(parallel)
  if(!is(DMRs, 'GRangesList')) DMRs <- split(DMRs, sign(DMRs$value))
  occs <- mclapply(DMRs, function(x) chromophobe::occupancy(y, x))
  return(occs)
} # }}}

## plotOccBySign -- plot occs and return them 
plotOccBySign <- function(occs=NULL, DMRs=NULL, y=NULL, dropQuiescent=F) { # {{{
  require(gridExtra)
  .grid.arrange.rows <- function(...) grid.arrange(..., nrow=1)
  if(is.null(occs)) occs <- getOccBySign(DMRs, y)
  changeNames <- c(gain='Gain of methylation', loss='Loss of methylation')
  if(dropQuiescent==TRUE) 
    changeNames <- paste(changeNames, '(exclusive of quiescent states)')
  names(occs) <- sub('^1', changeNames['gain'], 
                     sub('^-1', changeNames['loss'], 
                         names(occs)))
  plots <- lapply(names(occs), function(x) 
             plotOccupancy(occs[[x]],dropQuiescent=dropQuiescent)+ggtitle(x))
  do.call(.grid.arrange.rows, plots)
  return(occs)
} # }}}

