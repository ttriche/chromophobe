#' simple function for a donut plot of state width/count fractions 
#'
#' by default, use genomeFrac to determine how much of each state is represented
#'
#' @param  gs           a GenomicSegmentation or similar, or a data.frame
#' @param  text         text to put in the middle of the donut (NULL)
#' @param  colorScheme  named vector of colors (default is simplified REMC)
#' @param  ...          additional arguments to pass to genomeFrac
#' @param  fracSize     how big the text should be (5x typical label)
#' @param  legend       legend? (FALSE; no legend)
#'
#' @return a plot
#' 
#' @import ggplot2
#'
#' @export
donutPlot <- function(gs, text=NULL, colorScheme=NULL, ..., fracSize=5, legend=FALSE) { 

  if (is(gs, "data.frame")) {
    fracs <- gs
  } else {
    fracs <- genomeFrac(gs, ...)
  }
  if (is.null(colorScheme)) {
    colorScheme <- c("Promoter" = "red",
                     "Transcribed" = "darkgreen", 
                     "Enhancer" = "orange",
                     "Accessible" = "yellow",
                     "Het_Rpt_Qui" = "antiquewhite", 
                     "Bivalent" = "darkviolet",
                     "Repressed" = "gray50")
  }
  stopifnot(all(fracs$category %in% names(colorScheme)))
  
  p1 <- ggplot(fracs,
               aes(color=category,
                   fill=category, 
                   ymax=ymax, 
                   ymin=ymin, 
                   xmax=3, 
                   xmin=1)
               ) +
        geom_rect(color=NA, show.legend=legend) +
        scale_color_manual(values=colorScheme) +
        scale_fill_manual(values=colorScheme) +
        coord_polar(theta="y") +
        guides(color="none") +
        xlim(c(0, 3)) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.text= element_text(face="bold", family="serif", size=12),
          legend.title = element_text(face="bold", family="serif", size=16),
          panel.background = element_rect(fill=NA, color=NA),
        ) + 
        geom_text(
          aes(label=label, angle=angle, y=(ymin+ymax)/2, x=2), 
          inherit.aes=TRUE,
          color="white",
          fontface="bold",
          show.legend=FALSE,
          size=fracSize
        ) + 
        NULL
  
  if (!is.null(text)) {
    p1 <- p1 + annotate("text", 
                        x = 0, 
                        y = 0,
                        hjust = 0.5,
                        vjust = 0.5,
                        size = fracSize,
                        fontface = "bold", 
                        label = text)
  }

  if (legend) {
    p1 <- p1 + guides(fill=guide_legend(title=""))
  }

  return(p1)

}


# helper fn 
.addLabels <- function(dat, what=NULL, minpct=1) {
  
  dat$angle <- .getAngles(dat)
  dat$hjust <- .getHjust(dat$angle)
  dat$angle <- .fixAngles(dat$angle)
  dat$percent <- .getPercent(dat, minpct=minpct) 
  dat$label <- with(dat, .getLabel(percent, category, what=what))
  return(dat)

}


# helper fn 
.getAngles <- function(dat) {

  states <- nrow(dat)
  dat$end <- cumsum(dat$fraction) 
  dat$start <- c(0, dat$end[seq_len(states - 1)])
  dat$midpoint <- with(dat, ((start + end) / 2))
  return(90 - (360 * (dat$midpoint)))

}


# helper fn 
.fixAngles <- function(angle) ifelse( angle < -90, angle + 180, angle )


# helper fn
.getHjust <- function(angle) ifelse( angle < -90, 1, 0 )


# helper fn 
.getPercent <- function(dat, minpct=1) {

  pct <- round(dat$fraction * 100)
  pct <- ifelse(pct >= minpct, paste0(pct, "%"), "")
  return(pct)

}


# helper fn
.getLabel <- function(percent, category, what=NULL) {

  lab <- ifelse(percent != "", paste0(category, ": ", percent), "")
  if (!is.null(what)) lab <- paste(lab, "of" , what)
  return(lab)

}
