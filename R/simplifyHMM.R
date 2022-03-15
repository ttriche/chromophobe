#' convenience function for adding simplified names & colors to a ChromHMM track
#' 
#' @param HMM   a GRanges containing a ChromHMM track
#' @param cols  color table (will be loaded if not provided)
#' @param how   one of `MNEMONIC` (default), `STATE`, or `NUMBER`
#' 
#' @return      the same ChromHMM track, but simplified and colored simply
#'
#' @seealso     colorHMM
#'
#' @import      rtracklayer
#' 
#' @export
simplifyHMM <- function(HMM, cols=NULL, how=c("MNEMONIC","STATE","NUMBER")) { 

  if (is.null(cols)) {
    if (length(unique(HMM$name)) < 19) {
      message("Loading default Roadmap 18-state colors...")
      data(remc18state, package="chromophobe")
      cols <- remc18state
    } else { 
      message("Loading default Roadmap 25-state colors...")
      data(remc25state, package="chromophobe")
      cols <- remc25state
    }
  }
  
  if (how == "MNEMONIC" && !all(HMM$name %in% cols$MNEMONIC)) {
    # try stripping leading state numbers
    HMM$name <- sapply(strsplit(HMM$name, "_"), `[`, 2)
    stopifnot(all(HMM$name %in% cols$MNEMONIC))
  } else if (how == "NUMBER" && length(unique(HMM$name)) > nrow(cols)) {
    stop("More states in HMM than colors. Did you mean to use a different one?")
  }

  stopifnot(all(c("SIMPLE","RGBSIMPLE") %in% names(cols)))
  HMM$name <- cols[.matchState(HMM=HMM, cols=cols, how=how), "SIMPLE"]
  cols$RGB <- NULL
  names(cols)[which(names(cols) == "RGBSIMPLE")] <- "RGB"
  cols$MNEMONIC <- NULL
  names(cols)[which(names(cols) == "SIMPLE")] <- "MNEMONIC"
  return(colorHMM(HMM, cols=cols, how=how))

}
