bumps2bed <- function(bumps, pval=0.05, bedFile=NULL) { 
  require(bumphunter)
  require(rtracklayer)
  require(GenomicRanges)
  bumps.GR <- df2GR(bumps$table, keep=TRUE)
  bumps.GR <- bumps.GR[ bumps.GR$p.value <= pval ]
  bumps.GR$name<- paste0('bump', names(bumps.GR))
  bumps.GR$score <- bumps.GR$p.value
  export(bumps.GR, bedFile, "bed")
}
