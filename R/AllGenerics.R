## workhorse tabulation function
getOccupancy <- function(object, x=NULL) { # {{{
  if(is.null(x)) {
    states <- levels(mcols(object)[,'state'])
    byState <- split(object, mcols(object)[,'state'])
    unlist(lapply(byState, function(s) sum(as.numeric(width(ranges(s))))))
  } else { 
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
    getOccupancy(disjoint)
  }
} # }}}

setGeneric('occupancy', function(object,x,...) standardGeneric('occupancy'))
setMethod('occupancy', signature(object='JointSegmentation', x='missing'), #{{{
          function(object) {
            occ <- occupancy(rowData(object))
            states(occ) <- states(object)
            occ
          }) # }}}
setMethod('occupancy', signature(object='JointSegmentation', x='GenomicRanges'),#{{{
          function(object, x) {
            occ <- occupancy(rowData(object), x)
            states(occ) <- states(object)
            occ
          }) # }}}
setMethod('occupancy', signature(object='GRangesList', x='missing'),#{{{
          function(object) {
            states <- levels(mcols(object[[1]])[,'state'])
            d <- DataFrame(lapply(object, getOccupancy))
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

## the workhorses are these next two calls to getOccupancy()
setMethod('occupancy', signature(object='GRanges', x='missing'),#{{{
          function(object) getOccupancy(object)) # }}}
setMethod('occupancy', signature(object='GRanges', x='GRanges'),#{{{
          function(object, x) getOccupancy(object, x)) # }}}

## workhorse plot function
plotOccupancy <- function(x, dropQuiescent=FALSE, stateColors=NULL) { # {{{
  require(grid)
  require(ggplot2)
  title <- 'Chromatin state occupancy'
  x[,'state'] <- factor(rownames(x))
  if(dropQuiescent==TRUE) {
    x <- x[grep('quiescent', invert=TRUE, tolower(rownames(x))),]
    title <- paste(title, '(exclusive of quiescent states)')
  }
  if(!is.null(stateColors)) {
    message('automate this to handle states and colors from x')
    colorscale <- scale_color_manual(stateColors[unique(x[,'state'])])
  }
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
         ylab(paste0('Cumulative fraction of sites occupied')) + 
         xlab('') + 
         ggtitle(title)
  return(p)  
} # }}}

setMethod('plot', signature(x='Occupancy', y='missing'), # {{{
          function(x, ...) plotOccupancy(x)) # }}}
setMethod('plot', signature(x='JointSegmentation', y='GRanges'), # {{{
          function(x, y, ...) plotChromHMM(x, y, ...)) # }}}
setMethod('plot', signature(x='JointSegmentation', y='character'), # {{{
          function(x, y, ...) {
            if( y == 'emissions' ) {
              plotEmissions(emissions(x), ...)
            } else if( y == 'transitions' ) {
              plotTransitions(transitions(x), emissions(x), ...)
            } else if( y == 'occupancy' ) {
              plotChromHMM(x, ...)
            } else {
              stop(paste("Don't know how to plot",y,"for a JointSegmentation"))
            }
          }) # }}}
setMethod('plot', signature(x='JointSegmentation', y='missing'), # {{{
          function(x, ...) {
            message('Plotting occupancy; can also choose emissions/transitions')
            plot(occupancy(x), ...)
          }) # }}}

setGeneric('byState', function(object, ...) standardGeneric('byState'))
setMethod('byState', signature(object='Segmentation'), # {{{
          function(object) split(object, mcols(object)[,'state'])) # }}}

setGeneric('byChr', function(object, ...) standardGeneric('byChr'))
setMethod('byChr', signature(object='GenomicRanges'), # {{{ inherit from GR
          function(object) split(object, seqnames(object))) # }}}

setMethod('$', 'States', # {{{
  function(x, name) {
    xx <- x[[name, exact=FALSE]]
    names(xx) <- rownames(x)
    return(xx)
  }) # }}}

setGeneric('states', function(object, ...) { # {{{
           if(class(object) == 'States') return(object)
           return(slot(object, 'states'))
          }) # }}}

setGeneric('stateIds', function(object, ...) standardGeneric('stateIds'))
setMethod('stateIds', signature(object='ANY'), # {{{
          function(object) states(object)$Id) # }}}

setGeneric('stateNames', function(object, ...) standardGeneric('stateNames'))
setMethod('stateNames',  signature(object='ANY'), # {{{
          function(object) states(object)$Name) # }}}

setGeneric('stateColors', function(object, ...) standardGeneric('stateColors'))
setMethod('stateColors', signature(object='ANY'), # {{{
          function(object) states(object)$Color ) # }}}

toRgb <- function(color) { # {{{
  if(length(color) > 1) { 
    res <- t(sapply(color, toRgb))
    colnames(res) <- c('red','green','blue')
    return(res)
  } else {
    x <- try(col2rgb(color))
    if(!inherits(x, "try-error")) return(x)
    else return(NA)
  }
} # }}}
setGeneric('stateRGB', function(object, ...) standardGeneric('stateRGB'))
setMethod('stateRGB', signature(object='ANY'), # {{{
          function(object) toRgb(stateColors(object))) # }}}

toHex <- function(RGB) { # {{{
  if(is(RGB,'matrix')) {
    res <- apply(RGB, 1, toHex)
    names(res) <- rownames(RGB)
    return(res)
  } else if(is(RGB, 'character')) {
    return(toHex(toRgb(RGB)))
  } else if(length(RGB) == 3 && is.integer(RGB)) {
    return(rgb(red=RGB[1], green=RGB[2], blue=RGB[3], max=255))
  } else {
    return(NA)
  }
} # }}}
setGeneric('stateHex', function(object, ...) standardGeneric('stateHex'))
setMethod('stateHex', signature(object='ANY'), # {{{
          function(object) toHex(stateRGB(object))) # }}}

setGeneric('stateTypes', function(object,...) standardGeneric('stateTypes'))
setMethod('stateTypes',  signature(object='ANY'), # {{{
          function(object) states(object)$Type) # }}}

setGeneric('states<-', function(object, value, ...) # {{{
  standardGeneric('states<-'))  # }}}
setMethod('states<-', signature(object='JointSegmentation',value='States'), #{{{
          function(object, value) {
            message('This method ought to check states & reassign/re-name them')
            object@states <- value
            return(object)
          }) # }}}
setMethod('states<-', signature(object='Occupancy',value='States'), #{{{
          function(object, value) {
            if(all(!is.na(match(rownames(object), rownames(value))))) {
              object@states <- value
            } else if(all(!is.na(match(rownames(object), value$Name)))) {
              rownames(value) <- value[,'Name']
              object@states <- value
            } else if(all(!is.na(match(rownames(object), value$Id)))) {
              rownames(value) <- value[,'Id']
              object@states <- value
            } else {
              stop("Don't know how to map states without Names or Ids... !")
            }
            return(object)
          }) # }}}

setGeneric('probinit', function(object, ...) standardGeneric('probinit')) 
setMethod('probinit', signature(object='JointSegmentation'), # {{{
          function(object) object@probinit) # }}}

setGeneric('emissions', function(object, ...) standardGeneric('emissions'))
setMethod('emissions', signature(object='JointSegmentation'), # {{{
          function(object) object@emissions) # }}}

setGeneric('transitions', function(object, ...) standardGeneric('transitions'))
setMethod('transitions', signature(object='JointSegmentation'), # {{{
          function(object) object@transitions) # }}}

setGeneric('segmentation',function(object,x,...)standardGeneric('segmentation'))
setMethod('segmentation',signature(object='JointSegmentation',x='character'), # {{{
          function(object, x) rowData(object)[[x]]) # }}}
setMethod('segmentation',signature(object='JointSegmentation',x='missing'),#{{{
          function(object, x) rowData(object)) # }}}

setGeneric('posterior', function(object, x, ...) standardGeneric('posterior'))
setMethod('posterior', signature(object='JointSegmentation', x='character'),#{{{
          function(object, x) {
            stop('FIXME: Posterior probabilities are not yet supported!')
          }) # }}}

setMethod("combine", signature=signature(x="SummarizedExperiment", 
                                         y="SummarizedExperiment"), 
          function(x, y, ...) { # {{{
              if (class(x) != class(y)) {
                stop(paste("Error: objects must be the same class, but are ",
                           class(x), ", ", class(y), sep=""))
              }
              if( all(is.na(genome(rowData(x)))) || 
                  all(is.na(genome(rowData(y)))) ||
                  unique(na.omit(genome(rowData(x)))) !=
                  unique(na.omit(genome(rowData(y)))) ) {
                stop("Error: x and y have differing or unspecified genomes")
              }
              ## FIXME: allow for "packing out" missing features using NAs
              if( length(intersect(rownames(x), rownames(y))) < nrow(x) || 
                  length(intersect(rownames(x), rownames(y))) < nrow(y) ) {
                stop("Error: x and y have differing features, cannot combine")
              }
              commonAsys <- intersect(names(assays(x, withDimnames=FALSE)), 
                                      names(assays(y, withDimnames=FALSE)))
              if(length(commonAsys) < 1) stop('Error: no assays in common')
              else names(commonAsys) <- commonAsys

              combineAssay <- function(asy, x, y) {
                cbind(assays(x)[[asy]], assays(y[rownames(x), ])[[asy]])
              }

              SE <- SummarizedExperiment(
                assays=lapply(commonAsys, combineAssay, x=x, y=y),
                colData=merge(colData(x), colData(y), all=T, sort=F),
                rowData=rowData(x)
              )
              colnames(SE) <- c(colnames(x), colnames(y))
              return(SE)

          }) # }}}

setMethod("combine", signature=signature(x="JointSegmentation",
                                         y="JointSegmentation"), 
          function(x, y, ...) stop("Don't combine HMMs; reload them jointly."))

setMethod("keepSeqlevels", signature(x="SummarizedExperiment", value="ANY"),
          function(x, value) { # {{{
            y <- which(rownames(x) %in% names(keepSeqlevels(rowData(x),value)))
            return(x[y, ])
          }) # }}}

setMethod("sort", signature(x="SummarizedExperiment"),
          function(x) { # {{{ 
            x[ names(sort(rowData(x))), ]  
          }) # }}}

setGeneric("tabplot", function(dat, ...) standardGeneric('tabplot'))
## methods are defined in tabplot.R
