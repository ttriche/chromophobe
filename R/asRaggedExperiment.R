#' Take a set of segmentations and make them into a RaggedExperiment.
#' 
#' This seems pretty banal, and in some respects it is, but there are 
#' some fiddly bits to computing on segmentations, and this function 
#' attempts to ease them. For example, it simply isn't possible to use
#' sparseAssay or disjoinAssay with character, integer, or factor matrices.
#' If all states are represented in all segmentations, this reduces to 
#' generating and keeping a state-mapping key in metadata(), then referring
#' to it on the way out of the RaggedExperiment (usually this will employ 
#' disjoinSummarizedExperiment to map a bunch of regions across a bunch of 
#' segmentation models). Yet another possibility is that each segmentation 
#' has a probability for each possible state, and a collection of segmentations
#' would be better stored as a tensor (stack of arrays). This function attempts
#' to finesse the most painful parts of the user experience (based on mine).
#' Irritating problems like missing states can make life difficult, so this
#' function deals with them up-front by enumerating states, then keying them.
#' 
#' The function will probably be wrapped by a coercion before too much longer.
#' 
#' @param segs    segmentations, usually a GRangesList or GenomicSegmentations
#' @param scHMM   scChromHMM-like posterior probabilities? (different tricks)
#' @param keycol  mcols column to find the states that comprise a key ("name")
#' @param to      what to call this keycol to in the RaggedExperiment? ("state")
#' 
#' @return        a RaggedExperiment RE, with `key` and `yek` in metadata(RE)
#' 
#' @details       The `states` method uses `key` on RE assays.
#'                metadata(RE)$yek, as its name implies, reverses `key`.
#'
#' 
#' @examples 
#'
#' # some simplified chrY tracks 
#' data(simpleY, package="chromophobe")
#' asRaggedExperiment(simpleY)
#'
#' @seealso       states
#'
#' @import        RaggedExperiment
#' 
#' @export
asRaggedExperiment <- function(segs, scHMM=FALSE, keycol="name", to="state") {

  if (scHMM) stop("scChromHMM support is not yet complete.")
  if (elementType(segs) %in% c("GRanges", "GenomicSegmentation")) {
    key <- levels(factor(mcols(unlist(segs))[, keycol]))
  } else {
    key <- levels(factor(mcols(segs)[, keycol]))
  }
  yek <- seq_along(key)
  names(yek) <- key 
  
  RE <- RaggedExperiment(lapply(segs, .encode, yek=yek, keycol=keycol, to=to))  
  metadata(RE)$key <- key 
  metadata(RE)$yek <- yek 
  return(RE) 

}


# helper fn
.encode <- function(seg, yek, keycol="name", to="state") {

  mcols(seg)[, keycol] <- as.numeric(yek[mcols(seg)[, keycol]])
  mcols(seg) <- mcols(seg)[, keycol, drop=FALSE]
  names(mcols(seg)) <- to
  return(seg)

}


# helper fn for testing 
.decode <- function(asy, key) {

  res <- apply(asy, 2, function(x) key[x])
  rownames(res) <- rownames(asy)
  return(res)

}


# helper fn for testing
.recode <- function(asy, yek) {

  res <- apply(asy, 2, function(x) yek[x])
  rownames(res) <- rownames(asy)
  return(res)

}
