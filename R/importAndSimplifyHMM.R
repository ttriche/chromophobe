#' import and simplify a REMC or ENCODE ChromHMM track
#' 
#' @param x       a BED filename with a ChromHMM segmentation in it
#' @param genome  a genome associated with the segmentation (default is mm10)
#' @param species a species for the segmentation (default: guess from genome)
#' @param cols    a color key for simplifying ChromHMMs (pass to simplifyHMM)
#' @param verbose be verbose? (FALSE) 
#' @param ...     arguments (`only3`, `BPPARAM`) for aggregateStates
#' 
#' @return        a GRanges of aggregated, simplified HMM states
#' 
#' @import        rtracklayer
#' @import        GenomeInfoDb 
#' @import        GenomicRanges
#' 
#' @export
importAndSimplifyHMM <- function(x, genome="mm10", species=NULL, cols=NULL, verbose=FALSE, ...) {

  if (verbose) message("Importing and simplifying ", x, "... ", appendLF=FALSE)
  HMM <- rtracklayer::import(x)
  genome(HMM) <- genome
  seqlevelsStyle(HMM) <- "UCSC"
  mcols(HMM) <- DataFrame(name=mcols(HMM)[, "name"])
  if (is.null(species)) {
    species <- switch(genome,
                      mm9="Mus_musculus", 
                      mm10="Mus_musculus", 
                      GRCm38="Mus_musculus",
                      GRCm39="Mus_musculus",
                      hg18="Homo_sapiens", 
                      hg19="Homo_sapiens", 
                      hg38="Homo_sapiens",
                      GRCh37="Homo_sapiens",
                      GRCh38="Homo_sapiens",
                      dr11="Danio_rerio",
                      GRCz11="Danio_rerio")
  }
  HMM <- keepStandardChromosomes(HMM, species=species, pruning.mode="coarse")
  if (verbose) message("Done.")

  aggregateStates(simplifyHMM(HMM, cols=cols), verbose=verbose, ...)

}
