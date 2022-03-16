#' convenience function for generating donut plots
#' 
#' @param vec     vector of items
#' @param minpct  lower bound to label a slice
#'
#' @return a data.frame
#'
#' @import ggplot2
#' 
#' @export
makeFrac <- function(vec, minpct=7) {
  tbl <- table(vec, exclude=FALSE)
  dat <- data.frame(
    category=factor(names(tbl)),
    count=unname(as(tbl, "integer"))
  )
  dat$fraction <- dat$count / sum(dat$count)
  dat$percent <- round(dat$fraction * 100)
  dat$percent <- ifelse(dat$percent >= minpct, paste0(dat$percent, "%"), "")
  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n=-1)) 
  return(dat)
}
