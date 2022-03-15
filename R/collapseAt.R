#' Collapse a SummarizedExperiment-like object over a GRanges, sensibly.
#' Handy for things like chromophobe model testing, GenometriCorr, etc.
#' 
#' @param x         something descended from SummarizedExperiment, with assays()
#' @param y         something descended from a GRanges, with a granges() method
#' @param how       how to summarize the assay values (default: median)
#' @param BPPARAM   BPPARAM object to run in parallel (default: SerialParam())
#'
#' @return          an object with same colData but with new rowData and assays
#' 
#' @details         TODO: let the user provide their own summarizer function
#' 
#' @examples
#'
#' set.seed(123)
#' rr <- as(c("chr8:11657801-11665899","chr8:11665900-11700200"), "GRanges")
#' cd <- DataFrame(sample=letters[1:4], treatment=c(0,1,0,1))
#' rownames(cd) <- cd$sample
#' beta <- matrix(runif(8), nrow=length(rr), ncol=nrow(cd), 
#'                dimnames=list(as.character(rr), rownames(cd)))
#' tpm <- matrix(exp(rnorm(8)), nrow=length(rr), ncol=nrow(cd), 
#'               dimnames=list(as.character(rr), rownames(cd)))
#' se <- SummarizedExperiment(list(beta=beta,tpm=tpm), rowRanges=rr, colData=cd)
#' gr <- as("chr8:11665843-11665969", "GRanges")
#' collapseAt(se, gr)
#'
#' @import          GenomicRanges
#' @import          BiocParallel
#' @import          matrixStats
#' 
#' @export
collapseAt <- function(x, y, how=c('median','mean','sum'), BPPARAM=SerialParam()) {

  if (!identical(disjoin(y), y)) {
    message("You have overlapping target ranges. collapseAt will likely fail.")
  }
  how <- match.arg(how)
  fn <- switch(how,
               'median'=colMedians,
               'mean'=colMeans2,
               'sum'=colSums2)

  y <- subsetByOverlaps(granges(y), x)
  names(y) <- as.character(y) # coords
  x <- subsetByOverlaps(x, y)
  res <- x[seq_along(y), ] 
  ol <- findOverlaps(x, y)
  rowRanges(res) <- y

  # usually regions >> assays, so parallelize over regions rather than assays
  for (an in assayNames(x)) {
    assay(res, an) <- .coll(assay(x, an), ol=ol, y=y, fn=fn, BPPARAM=BPPARAM)
  }
  metadata(res)$how <- how
  return(res)

}


# helper fn
.coll <- function(asy, ol, y, fn, BPPARAM=SerialParam()) {

  byOverlap <- split(queryHits(ol), names(y)[subjectHits(ol)])
  res <- do.call(rbind, 
                 bplapply(byOverlap, 
                          function(i) fn(asy[i, , drop=F]), 
                          BPPARAM=BPPARAM))
  colnames(res) <- colnames(asy)
  rownames(res) <- names(y)
  return(res)

}
