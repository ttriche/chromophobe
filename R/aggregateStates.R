#' merge segments that are all in a row, but avoid slopping across chroms
#' 
#' @param HMM     a ChromHMM GRanges
#' @param only3   keep only `name`, `thick`, `itemRgb` as mcols? (FALSE)
#' @param verbose prattle endlessly about which contig we're on? (FALSE) 
#' @param BPPARAM BiocParallel object (default is SerialParam()) 
#' 
#' @return        HMM but without unnecessary segmentations (or mcols)
#' 
#' @import        GenomicRanges
#' @import        BiocParallel
#'
#' @export
aggregateStates <- function(HMM, only3=FALSE, verbose=FALSE, BPPARAM=SerialParam()) {

  grl <- split(HMM, seqnames(HMM))[unique(seqnames(HMM))]
  unname(unlist(GRangesList(bplapply(grl, 
                                     .Rle, 
                                     only3=only3, 
                                     verbose=verbose,
                                     BPPARAM=BPPARAM))))

}


# use an Rle to encode runs of states
# this function is imperfect but works OK
.Rle <- function(HMM, only3=FALSE, verbose=FALSE) { 

  chrom <- unique(seqnames(HMM)) 
  if (length(chrom) > 1) stop("This function should be run per-contig!")
  if (verbose) message("Aggregating states on ", chrom, "... ", appendLF=FALSE)
  columns <- c("name", "thick", "itemRgb")

  ungapped <- HMM
  mc <- names(mcols(HMM))
  mcd <- setdiff(columns, mc)
  mci <- intersect(columns, mc) 
  for (m in mcd) mcols(ungapped)[, m] <- NA
  names(ungapped) <- rep("HMM", length(ungapped))
  if (only3) mcols(ungapped) <- mcols(ungapped)[, columns]
  HHMM <- ungapped
  
  gapped <- gaps(HMM)
  ucd <- setdiff(mc, columns)
  if (length(gapped) > 0) {
    gapped$name <- "Gap" 
    gapped$thick <- ranges(gapped)
    gapped$itemRgb <- .toHex("0,0,0")
    names(gapped) <- rep("gaps", length(gapped))
    if (!only3) {
      for (u in ucd) {
        mcols(gapped)[, u] <- NA
      }
    }
    mcols(gapped) <- mcols(gapped)[, names(mcols(HHMM))]
    HHMM <- sort(c(ungapped, gapped))
  }

  asRle <- Rle(factor(HHMM$name))
  gr <- sort(HHMM[start(asRle)])
  end(gr) <- end(HHMM[end(asRle)])
  gr$name <- as.character(gr$name)
  gr$thick <- ranges(gr)

  gr <- unname(subset(gr, name != "Gap"))
  if (only3) {
    mcols(gr) <- mcols(gr)[, columns]
  }
  if (verbose) message("Done.")
  return(gr)

}
