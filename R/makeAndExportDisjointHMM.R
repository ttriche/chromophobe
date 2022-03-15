#' convenience function to disjoin and dump a set of HMMs
#' 
#' @param HMMs        a GRangesList of HMMs
#' @param trackName   what to name the track (will also generate BED filename)
#' @param BPPARAM     a SerialParam(), by default
#' 
#' @return            status of the export operation
#'
#' @import            rtracklayer
#' @import            BiocParallel
#' 
#' @export
makeAndExportDisjointHMM <- function(HMMs, trackName, BPPARAM=SerialParam()) {

  trackGenome <- unique(genome(HMMs))
  bedFile <- paste(gsub(" ", "_", trackName), trackGenome, "bed", sep=".")
  message("Writing disjoint HMM to ", bedFile, ".bgz ...")
  simpleCols <- .getSimpleCols()
  simpleCols$hex <- .toHex(simpleCols$RGB)
  dj <- .getDjHMM(HMMs, clarify=TRUE, BPPARAM=BPPARAM)
  dj$itemRgb <- simpleCols[dj$name, "hex"]
  mcols(dj) <- mcols(dj)[, c("name", "thick", "itemRgb")]
  dj$itemRgb[is.na(dj$itemRgb)] <- .toHex("250,250,250")
  trackLine <- new("TrackLine", name=trackName)
  export(dj, bedFile, trackLine=trackLine, index=TRUE)
  message("Exported and indexed ", bedFile, ".bgz")

}

