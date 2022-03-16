#' convenience function for generating donut plots
#' 
#' @param gs      a GRanges or GenomicSegmentation with states
#' @param minpct  lower bound (width-based) to label a state (1) 
#'
#' @return a data.frame
#'
#' @export
widthFrac <- function(gs, minpct=1) {

  states <- nlevels(factor(gs$name))
  bywidth <- vapply(split(gs, gs$name), function(x) sum(width(x)), integer(1))
  dat <- data.frame(
    category=factor(names(bywidth)),
    count=unname(bywidth)
  )
  dat$fraction <- dat$count / sum(dat$count)
  dat$percent <- round(dat$fraction * 100)
  dat$percent <- ifelse(dat$percent >= minpct, paste0(dat$percent, "%"), "")
  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n=-1)) 
  return(dat)

}

