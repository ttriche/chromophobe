testEnhancers <- function(HMM, p300=NULL, RNAseq=NULL, DHS=NULL, Pol2=NULL) {
  if(all(is.null(c(p300, RNAseq, DHS, Pol2)))){ 
    stop('You need to supply one or more of p300, RNAseq, DHS, or Pol2 data!')
  } else {
    require(GenometriCorr)
    stop('FIXME: Enhancer testing and regression is not yet finished')
  }
}
