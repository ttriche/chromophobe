#' simple function for a donut plot of segment width fractions 
#'
#' by default, use widthFrac to determine how much of each state is represented
#'
#' @param  gs           a GenomicSegmentation or similar 
#' @param  colorScheme  named vector of colors (default is simple REMC 25 state)
#' @param  ...          additional arguments to pass to widthFrac
#' @param  leg          legend? (FALSE; no legend)
#' @param  legTitle     title for legend? ("") 
#'
#' @return a plot
#' 
#' @import ggplot2
#'
#' @export
donutPlot <- function(gs, colorScheme=NULL, ..., leg=FALSE, legTitle="") {

  if (is.null(colorScheme)) {
    colorScheme <- c("Promoter" = "red",
                     "Transcribed" = "darkgreen", 
                     "Enhancer" = "orange",
                     "Accessible" = "yellow",
                     "Het_Rpt_Qui" = "white", 
                     "Bivalent" = "darkviolet",
                     "Repressed" = "gray")
  }
  fracs <- widthFrac(gs, ...)
  fracSize <- ifelse(min(fracs$fraction > 0.2), 9,
                     ifelse(min(fracs$fraction > 0.1), 8, 7))
  if (leg) fracSize <- 6 # accommodate legend
  p1 <- ggplot(fracs,
               aes(fill=category, ymax=ymax, ymin=ymin, xmax=3, xmin=1)) +
        geom_rect(color=NA, show.legend=leg) +
        coord_polar(theta="y") +
        xlim(c(0, 3)) +
        scale_fill_manual(values=colorScheme) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"),
          legend.text= element_text(face="bold", family="serif", size=12),
          legend.title = element_text(face="bold", family="serif", size=16),
          panel.background = element_rect(fill=NA, color=NA),
        )
  p1 <- p1 + geom_text(aes(label=percent, x=2, y=(ymin+ymax)/2),
                           color="white",
                           fontface="bold",
                           show.legend=FALSE,
                           inherit.aes=TRUE,
                           size=fracSize)
  if (leg) {
    p1 <- p1 + guides(fill=guide_legend(title=legTitle))
  }
  return(p1)

}
