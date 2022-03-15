# convenience function for rgb to hex 
.toHex <- function(rgbcsv) {
  
  rgbmat <- do.call(rbind, lapply(strsplit(rgbcsv, ","), as.integer))
  apply(rgbmat, 1, function(x) rgb(x[1], x[2], x[3], max=255))

}


# convenience function for adding colors to a ChromHMM track
addColors <- function(HMM, cols=NULL, how=c("MNEMONIC","STATE","NUMBER")) { 

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

  how <- match.arg(how)
  stopifnot(how == "MNEMONIC")
  rownames(cols) <- cols[, how]
  if (!all(HMM$name %in% rownames(cols))) {
    # try stripping leading state numbers
    HMM$name <- sapply(strsplit(HMM$name, "_"), `[`, 2)
    stopifnot(all(HMM$name %in% rownames(cols)))
  }
  
  # last minute conversion
  cols$hex <- .toHex(cols$RGB)
  HMM$itemRgb <- cols[HMM$name, "hex"]
  HMM$thick <- ranges(HMM)
  return(HMM)

}


# convenience function for adding simplified names & colors to a ChromHMM track
simplify <- function(HMM, cols=NULL, how=c("MNEMONIC","STATE","NUMBER")) { 

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
  stopifnot(all(c("SIMPLE","RGBSIMPLE") %in% names(cols)))
  cols$RGB <- cols$SIMPLERGB
  HMM$name <- cols$SIMPLE
  return(addColors(HMM, cols=cols, how="MNEMONIC"))

}
