#' simple function for a circular barchart of state transitions in a stackedHMM
#'
#' @param  stacked      a StackedHMM or similar (a GRanges of one is fine)
#' @param  text         text to put in the middle of the plot (NULL)
#' @param  colorScheme  named vector of colors (default is simplified REMC)
#' @param  legend       legend? (FALSE; no legend)
#'
#' @return a plot
#' 
#' @import ggplot2
#'
#' @export
stackedPlot <- function(stacked, text=NULL, colorScheme=NULL, legend=FALSE) { 

  if (is(stacked, "data.frame")) {
    fracs <- stacked  
  } else {
    fracs <- genomeFrac(stacked, ...) 
  }
  if (is.null(colorScheme)) {
    colorScheme <- c("Promoter" = "red",
                     "Transcribed" = "darkgreen", 
                     "Enhancer" = "orange",
                     "Accessible" = "yellow",
                     "Het_Rpt_Qui" = "antiquewhite", 
                     "Other_Unk" = "antiquewhite", 
                     "Bivalent" = "darkviolet",
                     "Repressed" = "gray50")
  }
  stopifnot(all(fracs$category %in% names(colorScheme)))
  
  fracSize <- ifelse(legend, 4, 5)
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
    p1 <- p1 + annotate(geom = 'text', x = 0.5, y = 0, label = text)
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
