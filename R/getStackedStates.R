#' pull out ranges from a stacked HMM that differ by condition (":" in names)
#' 
#' @param stacked   the stacked HMM 
#' @param state     the state of interest (default is "Bivalent") 
#' 
#' @return          a GRanges with the state for each condition in mcols() 
#' 
#' @details         Try subsetOverlaps(stacked, getStackedStates(stacked,state))
#'
#' @export
getStackedStates <- function(stacked, state="Bivalent") {

  hasState <- subset(stacked, grepl(paste0(state, ":"), name))
  pull <- subsetByOverlaps(stacked, hasState)
  mcols(pull) <- DataFrame(condition=sapply(strsplit(pull$name, ":"), `[`, 2),
                           state=sapply(strsplit(pull$name, ":"), `[`, 1))
  stopifnot(length(unique(pull)) == (length(pull)/2))
  res <- unique(granges(pull))
  for (cond in getConditions(hasState)) {
    mcols(res)[[cond]] <- subset(pull, condition == cond)$state
  }
  return(res)

}
