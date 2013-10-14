grSetToBsSeq <- function(grset, defaultCov=15) {
  require(bsseq)
  Cov <- matrix(rep(defaultCov, ncol(grset)*nrow(grset)), ncol=ncol(grset))
  M <- getBeta(grset) * defaultCov
  cnames <- colnames(grset)
  BSseq(M, Cov, pData=colData(grset), gr=granges(grset), sampleNames=cnames)
}

bedToBsSeq <- function(bedfile, MandCov='name', coefSmooth=NULL, genome='mm9'){
  require(bsseq)
  require(rtracklayer)
  grToBsSeq(import(bedfile, genome=genome))
}

grToBsSeq <- function(gr, MandCov='name', coefSmooth=NULL) {
  require(bsseq)
  message('This function assumes you have M/T as gr$',MandCov,' for each CpG.')
  MandT <- t(sapply(mcols(gr)[['name']], 
                    function(x) as.numeric(strsplit(x, '/')[[1]])))
  ## could be just as easily used for MethylSeekR import
  BSseq(M=MandT[,1,drop=FALSE], Cov=MandT[,2,drop=FALSE], gr=gr)
}

setAs('GenomicRatioSet','BSseq', function(from) grSetToBsSeq(from))
setAs('GenomicRanges','BSseq', function(from) grToBsSeq(from))
