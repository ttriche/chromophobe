setClassUnion('StatesOrNULL', c('NULL','States')) ## for loadChromHMM

## generic container for joint segmentations, inherits from SE
setClass('JointSegmentation', # {{{ a tweaked SummarizedExperiment
         representation(probinit='data.frame',
                        emissions='data.frame', 
                        transitions='matrix',
                        states='StatesOrNULL'),
         contains="SummarizedExperiment") # }}}

# setGeneric('states', function(object, ...) # {{{ (set in States-class.R) 
# standardGeneric('states'))  # }}}
setMethod('states', signature(object='JointSegmentation'), # {{{
          function(object) object@states) # }}}

setGeneric('states<-', function(object, value, ...) # {{{
  standardGeneric('states<-'))  # }}}
setMethod('states<-', signature(object='JointSegmentation',value='States'), #{{{
          function(object, value) {
            message('This method ought to check states & reassign/re-name them')
            object@states <- value
            return(object)
          }) # }}}

## iniial state probabilities
setGeneric('probinit', function(object, ...) standardGeneric('probinit')) 
setMethod('probinit', signature(object='JointSegmentation'), # {{{
          function(object) object@probinit) # }}}

## emission matrix
setGeneric('emissions', function(object, ...) standardGeneric('emissions'))
setMethod('emissions', signature(object='JointSegmentation'), # {{{
          function(object) object@emissions) # }}}

## transition matrix
setGeneric('transitions', function(object, ...) standardGeneric('transitions'))
setMethod('transitions', signature(object='JointSegmentation'), # {{{
          function(object) object@transitions) # }}}

## retrieve a hard-thresholded segmentation or subset thereof
setGeneric('segmentation',function(object,x,...)standardGeneric('segmentation'))
setMethod('segmentation',signature(object='JointSegmentation',x='character'), # {{{
          function(object, x) rowData(object)[[x]]) # }}}
setMethod('segmentation',signature(object='JointSegmentation',x='missing'),#{{{
          function(object, x) rowData(object)) # }}}

## posterior state probabilities: column per cell type, matrix per state?
setGeneric('posterior', function(object, x, ...) standardGeneric('posterior'))
setMethod('posterior', signature(object='JointSegmentation', x='character'),#{{{
          function(object, x) {
            stop('FIXME: Posterior probabilities are not yet supported!')
          }) # }}}

## genome-wide or GenomicRanges-wide occupancy for a JointSegmentation
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
