plotAsGviz <- function(object, GR, colors=NULL, ...) {
  require(Gviz)
  ## see https://stat.ethz.ch/pipermail/bioconductor/2013-January/050247.html
  ## e.g.
  ## Broad1 <- UcscTrack(track='Broad ChromHMM', 
  ##                     table="wgEncodeBroadHmmGm12878HMM", 
  ##                     trackType="AnnotationTrack", 
  ##                     genome='hg18', ## can be extracted from JointSeg object
  ##                     chromosome='chr18', 
  ##                     name='GM12878', 
  ##                     from=44675486, 
  ##                     to=44679944, 
  ##                     start="chromStart",  ## gratuitous for a GRanges
  ##                     end="chromEnd",      ## gratuitous for a GRanges
  ##                     feature="itemRgb", 
  ##                     id="name", 
  ##                     collapse=FALSE,
  ##                     stacking="dense")
  ## 
  ## feat <- unique(feature(Broad1))
  ## featCol <- setNames(as.list(rgb(t(sapply(strsplit(feat, ","), as.numeric)),
  ##                                 maxColorValue=255)), feat)
  ## displayPars(Broad1) <- featCol
  ## plotTracks(Broad1)
  ##
  stop('Plotting JointSegmentations in Gviz is almost, but not quite, done')
  plotTracks(JointSegAsTracks)
}
