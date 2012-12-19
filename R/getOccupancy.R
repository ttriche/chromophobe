getOccupancy <- function(object, x=NULL) {

  stopifnot(class(object) %in% c('Segmentation','GRanges'))

  if(is.null(x)) {
    unlist(lapply(split(object, mcols(object)[,'state']),
                  function(s) sum(as.numeric(width(ranges(s))))))
  } else { 
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
} 
