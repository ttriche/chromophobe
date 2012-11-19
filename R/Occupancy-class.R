setClass('Occupancy', contains="DataFrame")
Occupancy <- function(DF) { # {{{
  class(DF) <- 'Occupancy'
  return(DF)
} # }}}
setAs("Occupancy", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}

## relative occupancy for each state in a segmentation
setGeneric('occupancy', function(object,x,...) standardGeneric('occupancy'))
setMethod('occupancy', signature(object='JointSegmentation', x='missing'), #{{{
          function(object) occupancy(rowData(object))) # }}}
setMethod('occupancy', signature(object='JointSegmentation', x='GenomicRanges'),#{{{
          function(object, x) occupancy(rowData(object), x)) # }}}
setMethod('occupancy', signature(object='GRangesList', x='missing'),#{{{
          function(object) { 
            states <- levels(mcols(object[[1]])[,'state'])
#            stopifnot(all(unlist(lapply(object, function(x) 
#                          identical(states, levels(mcols(x)[,'state']))))))
            d <- DataFrame(lapply(names(object), function(m) 
                                  occupancy(object[[m]])))
            names(d) <- names(object)
            rownames(d) <- states
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ])
          }) # }}}
setMethod('occupancy', signature(object='GRangesList', x='GRanges'),#{{{
          function(object, x) { 
            states <- levels(mcols(object[[1]])[,'state'])
            stopifnot(all(unlist(lapply(object, function(x) 
                          identical(states, levels(mcols(x)[,'state']))))))
            d <- DataFrame(lapply(names(object), function(m) 
                                  occupancy(object[[m]], x)))
            names(d) <- names(object)
            rownames(d) <- states
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ])
          }) # }}}
setMethod('occupancy', signature(object='GRanges', x='missing'),#{{{
          function(object) { 
            stopifnot('state' %in% names(mcols(object)))
            states <- levels(mcols(object)[,'state'])
            n <- unlist(seqapply(split(object, mcols(object)$state), 
                                 function(s) sum(as.numeric(width(s)))))
            d <- DataFrame(occupancy=n/sum(n))
            rownames(d) <- states
            Occupancy(d)
          }) # }}}
setMethod('occupancy', signature(object='GRanges', x='GRanges'),#{{{
          function(object, x) { 
            stopifnot('state' %in% names(mcols(object)))
            ol <- findOverlaps(object, x)
            so <- object[queryHits(ol)]
            sx <- x[subjectHits(ol)]
            disjoint <- subsetByOverlaps(disjoin(c(sx,so,ignore.mcols=T)), so)
            mcols(disjoint)[,'state'] <- rep('', length(disjoint))
            byState <- split(so, mcols(so)[,'state'])
            for(state in names(byState)) {
              overlapping <- queryHits(findOverlaps(disjoint, byState[[state]]))
              if(length(overlapping)>0) {
                mcols(disjoint[overlapping])[,'state'] <- state
              }
            }
            mcols(disjoint)[,'state'] <- as.factor(mcols(disjoint)[,'state'])
            occupancy(disjoint)
          }) # }}}

## should speed things up... a little...
setMethod('plot', signature(x='Occupancy', y='missing'), # {{{
          function(x, ...) plotOccupancy(x)) # }}}

plotOccupancy <- function(x) { # {{{
  require(ggplot2)
  title <- paste('Occupancy by state')
  x[,'state'] <- factor(rownames(x))
  require(reshape2)
  byState <- melt(as.data.frame(x), id.vars='state')
  names(byState) <- gsub('^variable$', 'cell', names(byState))
  names(byState) <- gsub('^value$', 'fraction', names(byState))
  p <- ggplot(byState, aes(y=fraction, x=factor(cell), fill=state)) + 
         geom_bar(position="fill") + 
         theme(plot.title=element_text(face="bold",size=14),
               axis.title.x=element_text(face="bold",colour="#990000",size=14),
               axis.title.y=element_text(face="bold",colour="#990000",size=14),
               axis.text.x=element_text(angle=45, hjust=1, size=12), 
               axis.text.y=element_text(size=12)) +
         ylab(paste0('Cumulative fraction of sites occupied')) +
         xlab('Cell type') + 
         ggtitle(title)
  return(p)  
} # }}}
