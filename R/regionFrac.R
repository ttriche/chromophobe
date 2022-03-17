#' convenience function for generating donut plots of count proportions
#' 
#' @param gs      a GRanges or GenomicSegmentation with states
#' @param minpct  lower bound to label a slice
#' @param what    optional indicator
#'
#' @return a data.frame
#'
#' @export
regionFrac <- function(gs, minpct=1, what=NULL) {

  states <- nlevels(factor(gs$name))
  bycount <- vapply(split(gs, gs$name), function(x) length(x), integer(1))
  dat <- data.frame(category=factor(names(bycount)), count=unname(bycount))
 
  # stats 
  dat$fraction <- dat$count / sum(dat$count)
  dat$percent <- round(dat$fraction * 100)
  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n=-1)) 
 
  # add labels and angles for plotting 
  .addLabels(dat, what=what, minpct=minpct)

}
