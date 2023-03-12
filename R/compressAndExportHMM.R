#' Reduce a (usually simplified) segmentation, export to BED, compress, index.
#'
#' @param  HMM          a GenomicSegmentation or GRanges object
#' @param  name         a name to include as a trackLine in the BED file
#' @param  filename     the output (BED) file path (default: ./name.genome.bed)
#'
#' @import rtracklayer
#'
#' @examples 
#' 
#' data(chr19_HMM, package="chromophobe")
#' data(remc18state, package="chromophobe")
#' 
#' simpler <- remc18state
#' simpler$SIMPLE <- sub("(Promoter|Enhancer)", "Active", simpler$SIMPLE)
#' simpler[simpler$SIMPLE == "Active", "RGBSIMPLE"] <- "255,0,0" # Red
#' simpler$SIMPLE <- sub("(Transcribed|Het_Rpt_Qui)", "Other", simpler$SIMPLE)
#' simpler[simpler$SIMPLE == "Other", "RGBSIMPLE"] <- "255,255,255" # White
#' 
#' tmpdir <- tempdir()
#' setwd(tmpdir)
#' simplerHMM <- simplifyHMM(chr19_HMM, cols=simpler)
#' compressAndExportHMM(simplerHMM, name="chr19_HMM_simple", 
#'                      filename=file.path(tmpdir, "chr19_HMM.simple.mm10.bed"))
#' list.files(tmpdir) 
#'
#' @export
compressAndExportHMM <- function(HMM, name, filename=NULL) { 
  
  trackGenome <- unique(genome(HMM))
  trackLine <- new("TrackLine", name=name)
  if (is.null(filename)) filename <- paste(name, trackGenome, "bed", sep=".")
  export(aggregateStates(HMM), filename, index=TRUE) 

}
