#' list all the states found in a BED file of an HMM, stacked or otherwise 
#' 
#' @param HMM     the HMM 
#' @param asTable return a table of states? (FALSE) 
#' 
#' @return        depending on the value of asTable, either a table or a list
#' 
#' @export
getStates <- function(HMM, asTable=FALSE) { 

  tbl <- table(sapply(strsplit(HMM$name, ":"), `[`, 1))
  if (asTable) return(tbl)
  return(names(tbl))

}
