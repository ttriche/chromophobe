addHistoneTracks <- function(segs, 
                             path=NULL,
                             sep='-',
                             pattern='.pval.signal.bw',
                             marks=c('H3K27me3',
                                     'H3K36me3',
                                     'H3K4me1',
                                     'H3K27ac',
                                     'H3K4me3'), ...) {
  require('Gviz')
  addBigWigs <- function(segs) {
    for(track in names(segs)) {
      for(mark in marks) {
        bw <- paste0(mark,'.BigWig')
        colData(segs)[track, bw] <- paste0(track, sep, mark, pattern)
        if(!is.null(path)) {
          colData(segs)[track, bw] <- paste0(path, '/', colData(segs)[track, bw])
        }
      }
    }
  }
}
