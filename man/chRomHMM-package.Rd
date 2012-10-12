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
  model how many changes fall into a specific state it is a bit of a pain.
  This package makes it easier to relabel, prune/collapse, plot, and annotate
  the output of a ChromHMM run, runs, or other experiments mapped onto HMMs. 
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
  TY  - JOUR
  AU  - Ernst, Jason
  AU  - Kellis, Manolis
  TI  - ChromHMM: automating chromatin-state discovery and characterization
  JA  - Nat Meth
  PY  - 2012/03//print
  VL  - 9
  IS  - 3
  SP  - 215
  EP  - 216
  PB  - Nature Publishing Group, a division of Macmillan Publishers Limited. All Rights Reserved.
  SN  - 1548-7091
  UR  - http://dx.doi.org/10.1038/nmeth.1906
  M3  - 10.1038/nmeth.1906
  N1  - 10.1038/nmeth.1906
  L3  - http://www.nature.com/nmeth/journal/v9/n3/abs/nmeth.1906.html#supplementary-information
}
\keyword{ package }
\examples{

## not run:
## 
## import.ChromHMM(
##   file="~/Dropbox/ChromHMM/RoadCODE3.model25/BLADDER_25_segments.bed", 
##   as='GRL', 
##   genome='hg19')

}
