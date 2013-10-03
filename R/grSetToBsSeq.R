grSetToBsSeq <- function(grset, defaultCov=15) {
  require(bsseq)
  Cov <- matrix(rep(defaultCov, ncol(grset)*nrow(grset)), ncol=ncol(grset))
  M <- getBeta(grset) * defaultCov
  cnames <- colnames(grset)
  BSseq(M, Cov, pData=colData(grset), gr=granges(grset), sampleNames=cnames)
}

setAs('GenomicRatioSet','BSseq', function(from) grSetToBsSeq(from))
