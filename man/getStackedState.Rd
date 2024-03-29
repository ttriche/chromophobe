% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getStackedState.R
\name{getStackedState}
\alias{getStackedState}
\title{pull out ranges from a stacked HMM that differ by condition (":" in names)}
\usage{
getStackedState(stacked, state = "Bivalent")
}
\arguments{
\item{stacked}{the stacked HMM}

\item{state}{the state of interest (default is "Bivalent")}
}
\value{
a GRanges with the state for each condition in mcols()
}
\description{
pull out ranges from a stacked HMM that differ by condition (":" in names)
}
\details{
Try subsetOverlaps(stacked, getStackedState(stacked,state))
}
