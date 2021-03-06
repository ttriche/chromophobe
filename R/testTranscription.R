testTranscription <- function(HMM, CAGE=NULL, RNAseq=NULL, DHS=NULL, Pol2=NULL){
  if(all(is.null(c(RNAseq, CAGE, DHS, Pol2)))){ 
    stop('You need to supply one or more of CAGE, RNAseq, DHS, or Pol2 data!')
  } else {
    require(GenometriCorr)
    stop('FIXME: Transcription prediction testing is not yet finished')
  }
}
