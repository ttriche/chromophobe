% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addRoadmapColors.R
\name{addRoadmapColors}
\alias{addRoadmapColors}
\title{convenience function to relabel states with standard colors}
\usage{
addRoadmapColors(HMM, what = c("mnemonic", "state", "num"), cols = NULL)
}
\arguments{
\item{HMM}{a GRanges with !is.null(mcols(HMM)$name)}

\item{what}{is $name a "state", state"num", or "mnemonic"? ("mnemonic")}

\item{cols}{which colors to use? (default is 18-state Roadmap colors)}
}
\value{
a GRanges with adjusted colors
}
\description{
convenience function to relabel states with standard colors
}
