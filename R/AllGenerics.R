## Occupancy methods
getOccupancy <- function(object, x=NULL) { # {{{

  stopifnot(class(object) %in% c('Segmentation','GRanges'))

  if(is.null(x)) {
    ## FIXME: plot occupancy for a single segmentation needs work
    unlist(lapply(split(object, mcols(object)[,'state']),
                  function(s) sum(as.numeric(width(ranges(s))))))
  } else { 
    ## FIXME: there has to be a much faster way to do this
    if(!class(object) %in% c('Segmentation','GRanges')) occupancy(object, x)
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
            occ <- occupancy(segmentation(object))
            occ@states <- states(object)
            occ
          }) # }}}
setMethod('occupancy', signature(object='JointSegmentation', x='GenomicRanges'),#{{{
          function(object, x) {
            occ <- occupancy(segmentation(object), x)
            occ@states <- states(object)
            occ
          }) # }}}
setMethod('occupancy', signature(object='JointSegmentation', x='SummarizedExperiment'),#{{{
          function(object, x) {
            occ <- occupancy(segmentation(object), rowData(x))
            occ@states <- states(object)
            occ
          }) # }}}
setMethod('occupancy', signature(object='SegmentationList', x='missing'),#{{{
          function(object) {
            l <- lapply(object, occupancy)
            s <- unique(do.call(c, lapply(l, names)))
            d <- DataFrame(lapply(l, function(ll) {
              ll[setdiff(s, names(ll))] <- 0
              return(ll[s])
            }))
            rownames(d) <- s
            sts <- states(object)
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ], sts)
          }) # }}}
setMethod('occupancy', signature(object='SegmentationList', x='GRanges'),#{{{
          function(object, x) { 
            l <- lapply(object, occupancy, x=x)
            s <- unique(do.call(c, lapply(l, names)))
            d <- DataFrame(lapply(l, function(ll) {
              ll[setdiff(s, names(ll))] <- 0
              return(ll[s])
            }))
            rownames(d) <- s
            sts <- states(object)
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ], sts)
          }) # }}}
setMethod('occupancy', signature(object='GRangesList', x='missing'),#{{{
          function(object) {
            l <- lapply(object, occupancy)
            s <- unique(do.call(c, lapply(l, names)))
            d <- DataFrame(lapply(l, function(ll) {
              ll[setdiff(s, names(ll))] <- 0
              return(ll[s])
            }))
            rownames(d) <- s
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ])
          }) # }}}
setMethod('occupancy', signature(object='GRangesList', x='GRanges'),#{{{
          function(object, x) { 
            l <- lapply(object, occupancy, x=x)
            s <- unique(do.call(c, lapply(l, names)))
            d <- DataFrame(lapply(l, function(ll) {
              ll[setdiff(s, names(ll))] <- 0
              return(ll[s])
            }))
            rownames(d) <- s
            Occupancy(d[ order(rowMeans(as.matrix(d)), decreasing=TRUE), ])
          }) # }}}
setMethod('occupancy', signature(object='Segmentation', x='missing'),#{{{
          function(object) {
            occ <- getOccupancy(object)
            d <- DataFrame(occ)
            rownames(d) <- names(occ)
            Occupancy(d[ order(occ, decreasing=TRUE), ], states=states(object))
          }) # }}}
setMethod('occupancy', signature(object='Segmentation', x='GRanges'),#{{{
          function(object, x) {
            occ <- getOccupancy(object, x)
            d <- DataFrame(occ)
            rownames(d) <- names(occ)
            Occupancy(d[ order(occ, decreasing=TRUE), ], states=states(object))
          }) # }}}
setMethod('occupancy', signature(object='GRanges', x='missing'),#{{{
          function(object) getOccupancy(object)) # }}}
setMethod('occupancy', signature(object='GRanges', x='GRanges'),#{{{
          function(object, x) getOccupancy(object, x)) # }}}
setMethod('plot', signature(x='Occupancy', y='missing'), # {{{
          function(x, ...) plotOccupancy(x)) # }}}

## utility splitting methods
setGeneric('byState', function(object, ...) standardGeneric('byState'))
setMethod('byState', signature(object='Segmentation'), # {{{
          function(object) split(object, mcols(object)[,'state'])) # }}}
setGeneric('byChr', function(object, ...) standardGeneric('byChr'))
setMethod('byChr', signature(object='GenomicRanges'), # {{{
          function(object) split(object, seqnames(object))) # }}}

## States methods
setGeneric('states', function(object, ...) { # {{{
           if(class(object) == 'States') {
             return(object)
           } else if('states' %in% slotNames(object)) {
             return(slot(object, 'states'))
           } else {
             return(NULL)
           }
          }) # }}}
setMethod('$', 'States', # {{{
  function(x, name) {
    xx <- x[[name, exact=FALSE]]
    names(xx) <- rownames(x)
    return(xx)
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
setMethod('stateColors', signature(object='Occupancy'), # {{{
          function(object) {
            colrs <- states(object)$Color 
            names(colrs) <- states(object)$Name
            return(colrs[rownames(object)])
          }) # }}}
setGeneric('stateTypes', function(object,...) standardGeneric('stateTypes'))
setMethod('stateTypes',  signature(object='ANY'), # {{{
          function(object) states(object)$Type) # }}}
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
setGeneric('states<-', function(object, value, ...) # {{{
  standardGeneric('states<-'))  # }}}
setMethod('states<-', signature(object='JointSegmentation',value='States'),# {{{
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

## JointSegmentation methods 
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
setMethod('segmentation',signature(object='JointSegmentation',x='missing'),#{{{
           function(object, x) rowData(object)
          ) # }}}
setGeneric('posterior', function(object, x, ...) standardGeneric('posterior'))
setMethod('posterior', signature(object='JointSegmentation', x='character'),#{{{
          function(object, x) {
            stop('FIXME: Posterior probabilities are not yet supported!')
          }) # }}}
setMethod("combine", signature=signature(x="SummarizedExperiment", 
                                         y="SummarizedExperiment"),  # {{{
          function(x, y, ...) {
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
                                         y="JointSegmentation"), # {{{ 
          function(x,y,...) stop("Don't combine HMMs; reload them jointly"))#}}}
setMethod("keepSeqlevels", signature(x="SummarizedExperiment", value="ANY"),
          function(x, value) { # {{{
            y <- which(rownames(x) %in% names(keepSeqlevels(rowData(x),value)))
            return(x[y, ])
          }) # }}}
setMethod("sort", signature(x="SummarizedExperiment"),
          function(x) { # {{{ 
            x[ names(sort(rowData(x))), ]  
          }) # }}}
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

## SegmentationList methods (mostly to produce a valid Segmentation with states)
setMethod('$', 'SegmentationList', # {{{
  function(x, name) Segmentation(x[[name]], states(x))) # }}}

## methods defined in tabplot.R
setGeneric("tabplot", function(dat, ...) standardGeneric('tabplot'))
