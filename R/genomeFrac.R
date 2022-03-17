#' convenience function for generating donut plots of basepair proportions
#' 
#' @param gs      a GRanges or GenomicSegmentation with states
#' @param minpct  lower bound (width-based) to label a state (1) 
#' @param what    optional indicator
#'
#' @return a data.frame
#'
#' @export
genomeFrac <- function(gs, minpct=1, what=NULL) {

  states <- nlevels(factor(gs$name))
  bywidth <- vapply(split(gs, gs$name), function(x) sum(width(x)), integer(1))
  dat <- data.frame(category=factor(names(bywidth)), count=unname(bywidth))

  # stats
  dat$fraction <- dat$count / sum(dat$count)
  dat$percent <- round(dat$fraction * 100)
  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n=-1)) 

  # add labels and angles for plotting 
  .addLabels(dat, what=what, minpct=minpct)

}
