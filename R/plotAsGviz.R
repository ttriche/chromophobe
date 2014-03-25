asGvizTracks <- function(HMM, ranges, ...) {
  cols <- as.list(stateColors(HMM))
  tracks <- lapply(segmentations(HMM), 
                   function(x) {
                     xx <- subsetByOverlaps(x, ranges)
                     track <- AnnotationTrack(xx,
                                              stacking="dense",
                                              from=start(ranges),
                                              to=end(ranges),
                                              id=xx$state,
                                              feature=xx$state,
                                              background.title='orange')
                     displayPars(track) <- cols
                     return(track)
                   })
  for(i in names(tracks)) tracks[[i]]@name <- i
  return(tracks)
}

plotAsGviz <- function(object, ranges, colors=NULL, ...) {
  require(Gviz)
  plotTracks(asGvizTracks(object, ranges, ...))
}
