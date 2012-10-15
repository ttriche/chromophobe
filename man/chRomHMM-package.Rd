\name{chRomHMM-package}
\alias{chRomHMM-package}
\alias{chRomHMM}
\docType{package}
\title{
Package to plot genomic enrichments, as with ChromHMM models
}
\description{
  ChromHMM is a wonderful tool for heavy lifting and modeling chromatin states.
  As a post-processing tool it is slightly less wonderful, and if you want to 
  plot how cell-specific changes fall into various states, it's a bit of a pain.
  This package makes it easier to relabel, prune/collapse, plot, and/or annotate
  the output of a ChromHMM run, runs, or other experiments mapped onto segments.
}
\details{
\tabular{ll}{
Package: \tab chRomHMM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-10-12\cr
License: \tab APL\cr
}
More here soon.  This is a 0.1 release.  Try import.ChromHMM(yourBedFile).

}
\author{
Tim Triche, Jr. <tim.triche@gmail.com>
}
\references{
  ChromHMM: automating chromatin-state discovery and characterization 
  Jason Ernst & Manolis Kellis \cr
}
\keyword{ package }
\examples{

## not run:
## 
## bladder.HMM <- import.ChromHMM(
##   file="~/Dropbox/ChromHMM/RoadCODE3.model25/BLADDER_25_segments.bed", 
##   as='GRL', 
##   genome='hg19'
## )
##
## bladder <- import.posterior(
##   HMM=bladder.HMM, 
##   path="~/Dropbox/ChromHMM/RoadCODE3.model25/POSTERIOR",
##   states="~/Dropbox/ChromHMM/RoadCODE3.model25/labelStates.txt"
## )
## 
## plot.ChromHMM( 
##   bladder, 
##   what=c('genomic')
## )
## 

}
