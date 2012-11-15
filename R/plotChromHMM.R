## plot genome or site % occupancy for a given HMM's states
##
## FIXME: handle WIG tracks and other things with scores/depths, like segtools!
##
plotChromHMM <- function(object,GR=NULL,colors=NULL,asReads=F,name=NULL,...) {
  if(is.null(GR)) plot(occupancy(object))
  else plot(occupancy(object, GR))
}
