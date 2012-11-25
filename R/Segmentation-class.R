setClass('Segmentation', contains="GRanges")

setAs("GRanges", "Segmentation", function(from) { # {{{
  class(from) <- 'Segmentation'
  return(from)
}) # }}}

setGeneric('byState', function(object, ...) standardGeneric('byState'))
setMethod('byState', signature(object='Segmentation'), # {{{
          function(object) split(object, mcols(object)[,'state'])) # }}}

setGeneric('byChr', function(object, ...) standardGeneric('byChr'))
setMethod('byChr', signature(object='GenomicRanges'), # {{{ inherit from GR
          function(object) split(object, seqnames(object))) # }}}
