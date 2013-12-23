bumps2bed <- function(bumps, bedFile=NULL) { 
  require(bumphunter)
  require(rtracklayer)
  require(GenomicRanges)
  bumps.GR <- df2GR(bumps$table, keep=TRUE)
  bumps.GR$name<- paste0('bump', names(bumps.GR))
  bumps.GR$score <- bumps.GR$p.value
  export(bumps.GR, bedFile, "bed")
}
