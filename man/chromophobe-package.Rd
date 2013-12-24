\name{chromophobe-package}
\alias{chromophobe-package}
\alias{chromophobe}
\docType{package}
\title{
  Package to explore genomic enrichments for ChromHMM and related state models
}
\description{
  ChromHMM is a wonderful tool for heavy lifting and modeling chromatin states.
  As a post-processing tool it is slightly less wonderful, and if to plot how 
  cell-type-specific changes fall into various states, it's a bit of a pain.
  This package makes it easier to relabel, prune/collapse, plot, and/or annotate
  the output of a ChromHMM run, runs, or other experiments mapped onto segments.
  Functions are also provided to ease compositional analysis of DNA methylation
  data (see ?composition), import/export segmentations (see ?bumps2bed), and 
  various other conveniences often desired when working with these models. 
}
\details{
\tabular{ll}{
Package: \tab chromophobe\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-12-22\cr
License: \tab APL\cr
}
More here soon.  This is a 0.1 release.  Try import.ChromHMM(yourBedFile).
See also ?composition and ?bumps2bed for some convenience functions.
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
