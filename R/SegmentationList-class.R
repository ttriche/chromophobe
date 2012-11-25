setClass("SegmentationList", # {{{ a GRangesList w/Segmentation elements 
    prototype(elementType="Segmentation"),
    contains="GRangesList"
) # }}}

setAs("GRangesList", "SegmentationList", function(from) { # {{{
  from <- seqapply(from, function(x) {
    class(x) <- 'Segmentation'
    return(x)    
  }) 
  class(from) <- 'SegmentationList'
  return(from)
}) # }}}

SegmentationList <- function(...) { # {{{
    as(GRangesList(list(...)), 'SegmentationList')
} # }}}
