asGvizTracks <- function(HMM, ranges, ...) {
  cols <- as.list(stateColors(HMM))
  tracks <- lapply(segmentations(HMM), 
                   function(x) {
                     xx <- sort(subsetByOverlaps(x, ranges))
                     ## truncate so that all segments get plotted properly
                     start(xx[1]) <- start(sort(ranges)[1])
                     end(xx[length(xx)]) <- end(sort(ranges)[length(ranges)])
                     track <- AnnotationTrack(xx,
                                              stacking="dense",
                                              from=start(ranges),
                                              to=end(ranges),
                                              id=xx$state,
                                              feature=xx$state,
                                              background.title='darkred')
                     displayPars(track) <- append(cols, list(lwd.border=0))
                     return(track)
                   })
  for(i in names(tracks)) tracks[[i]]@name <- paste0(i, ' HMM')
  return(tracks)
}

plotAsGviz <- function(object, ranges, colors=NULL, ...) {
  require(Gviz)
  plotTracks(asGvizTracks(object, ranges, ...))
}
