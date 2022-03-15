#' convenience function to relabel states with standard colors 
#' 
#' @param HMM   a GRanges with !is.null(mcols(HMM)$name)
#' @param what  is $name a "state", state"num", or "mnemonic"? ("mnemonic")
#' @param cols  which colors to use? (default is 18-state Roadmap colors)
#' 
#' @return      a GRanges with adjusted colors
#' 
#' @export
addRoadmapColors <- function(HMM, what=c("mnemonic","state","num"), cols=NULL) {
 
  if (is.null(cols)) {
    message("Defaulting to 18-state Roadmap colors...")
    data("remc18state", package="chromophobe")
    (cols <- remc18state)
  }

  stop("Not finished")

}
