#' list all the conditions found in a BED file of a stacked HMM
#' 
#' @param stacked the stacked HMM 
#' 
#' @return        the conditions found in the stacked HMM (via state:condition)
#' 
#' @export
getConditions <- function(stacked) {

  levels(factor(sapply(strsplit(stacked$name, ":"), `[`, 2)))

}
