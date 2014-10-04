plotHistoneTracks <- function(jointSeg, 
                              track, 
                              trackname='ChIPseq', 
                              genome='hg19',
                              marks=c('H3K27me3', 'H3K36me3', 'H3K4me1',
                                      'H3K27ac', 'H3K4me3'), 
                              markColors=c(H3K27me3='gray50', 
                                           H3K36me3='green',
                                           H3K4me1='orange', 
                                           H3K27ac='darkorange',
                                           H3K4me3='darkred'),
                              ...) 
{

  splt <- function(x, y='/') strsplit(x, y, fixed=T)[[1]]

  strcut <- function(x, y='.', z=1) splt(x, y)[z]

  strpop <- function(x, y='/') splt(x, y)[length(splt(x, y))]

  getMark <- function(x) ifelse(length(x) > 1, # recurse if vector
                                  sapply(x, getMark), # otherwise,
                                  strcut(strcut(x, '-', 2), '.'))

  getColor <- function(wigname) markColors[getMark(wigname)]

  bigWigCols <- paste(marks, 'BigWig', sep='.')
  stopifnot(all(bigWigCols %in% colnames(colData(jointSeg))))
  names(bigWigCols) <- marks

  require(Gviz) ## since it is not yet a Depends:
  bigWigs <- colData(jointSeg)[track, bigWigCols]
  names(bigWigs) <- names(bigWigCols)
  tracks <- lapply(bigWigs, 
                   function(z) {
                     DataTrack(range=z,
                               size=1,
                               alpha=0.75,
                               type=c('h'),
                               genome=genome, 
                               name=trackname, 
                               background.title='saddlebrown',
                               fill=getColor(strpop(strcut(z))),
                               col=getColor(strpop(strcut(z))))
                   })
  OverlayTrack(tracks, background.title='saddlebrown', name=trackname)

}
