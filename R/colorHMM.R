#' add (new) colors to a ChromHMM track
#'
#' @param HMM   a GRanges from an ChromHMM segmentation
#' @param cols  a data.frame of colors and associated states
#' @param how   how to match the HMM states to those in cols? ("MNEMONIC")
#' 
#' @return      same HMM but with better colors
#' 
#' @import      GenomicRanges
#' @import      rtracklayer
#'
#' @export 
colorHMM <- function(HMM, cols=NULL, how=c("MNEMONIC","STATE","NUMBER")) { 

  if (is.null(cols)) {
    if (length(unique(HMM$name)) < 13) {
      message("Loading Blueprint 12-state colors...")
      data(blueprint12state, package="chromophobe")
      cols <- blueprint12state 
    } else if (length(unique(HMM$name)) < 19) {
      message("Loading default Roadmap 18-state colors...")
      data(remc18state, package="chromophobe")
      cols <- remc18state
    } else { 
      message("Loading default Roadmap 25-state colors...")
      data(remc25state, package="chromophobe")
      cols <- remc25state
    }
  }
  stopifnot("RGB" %in% names(cols))

  how <- match.arg(how)
  if (how == "MNEMONIC" && !all(HMM$name %in% cols$MNEMONIC)) {
    # try stripping leading state numbers
    HMM$name <- sapply(strsplit(HMM$name, "_"), `[`, 2)
    stopifnot(all(HMM$name %in% cols$MNEMONIC))
  } else if (how == "NUMBER" && length(unique(HMM$name)) > nrow(cols)) {
    stop("More states in HMM than colors. Did you mean to use a different one?")
  }
  cols$hex <- .toHex(as.character(cols$RGB))
  HMM$itemRgb <- cols[.matchState(HMM=HMM, cols=cols, how=how),  "hex"]
  HMM$thick <- ranges(HMM)
  return(HMM)

}


# convenience function for rgb to hex 
.toHex <- function(rgbcsv) {
  
  rgbmat <- do.call(rbind, lapply(strsplit(rgbcsv, ","), as.integer))
  apply(rgbmat, 1, function(x) rgb(x[1], x[2], x[3], max=255))

}


# convenience function for state matching 
.matchState <- function(HMM, cols, how=c("MNEMONIC","STATE","NUMBER")) {

  how <- match.arg(how)
  if (how == "MNEMONIC") idx <- match(HMM$name, cols$MNEMONIC)
  if (how == "STATE") idx <- match(HMM$name, rownames(cols))
  if (how == "NUMBER") idx <- as.integer(HMM$name)
  return(idx) 

}
