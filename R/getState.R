#' pull out ranges from a GenomicSegmentation or similar
#' 
#' trivial function but still handy 
#'
#' @param gs        the GenomicSegmentation or GRanges or similar
#' @param state     the state of interest (default is "Bivalent") 
#' 
#' @return          that subset of the GenomicSegmentation (or similar) 
#' 
#' @export
getState <- function(gs, state="Bivalent") {

  subset(gs, grepl(state, gs$name))

}
