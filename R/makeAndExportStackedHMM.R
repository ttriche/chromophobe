#' convenience function for dumping sets of chromHMM "contrasts" as stacked HMMs
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
makeAndExportStackedHMM <- function(HMMs, trackName, BPPARAM=SerialParam()) {

  trackGenome <- unique(genome(HMMs))
  bedFile <- paste(gsub(" ", "_", trackName), trackGenome, "bed", sep=".")
  message("Writing stacked HMM to ", bedFile, ".bgz ...")
  simpleCols <- .getSimpleCols()
  simpleCols$hex <- .toHex(simpleCols$RGB)
  stacked <- stackHMMs(HMMs, BPPARAM=BPPARAM)
  mcols(stacked) <- mcols(stacked)[, c("name", "thick", "itemRgb")]
  trackLine <- new("TrackLine", name=trackName)
  export(stacked, bedFile, trackLine=trackLine, index=TRUE)
  message("Exported and indexed ", bedFile, ".bgz")

}
