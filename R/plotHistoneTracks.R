plotHistoneTracks <- function(jointSeg, track, trackname='ChIPseq', genome='hg19',
                              marks=c('H3K27me3', 'H3K36me3', 'H3K4me1',
                                      'H3K27ac', 'H3K4me3') , ...) 
{

  strcut <- function(x, y='.', z=1) strsplit(x, y, fixed=T)[[1]][z]
  strpop <- function(x, y='/') {
    splt <- strsplit(x, y, fixed=T)[[1]]
    splt[length(splt)] 
  }
  getMark <- function(x) {
    if(length(x) > 1) sapply(x, getMark) 
    else strcut(strcut(x, '-', 2), '.')
  }
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
