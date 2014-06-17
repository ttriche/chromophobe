\name{bumps2bed}
\alias{bumps2bed}
\title{
  Grab output from bumphunter, filter by pval (default 0.05), export for Homer. 
}
\description{
  The homer package is perhaps the easiest way to assess motif enrichment for 
  defined intervals such as chromatin states or differentially enriched regions
  (as with DMRs, DiffBind regions, blocks, variable methylation regions, etc.).
  This function takes the output of bumphunter or blockFinder and dumps it into
  a BED file with columns suitably formatted so that Homer will not puke on it.
}
\details{
  See also ?composition for more convenience functions.
}
\author{
  Tim Triche, Jr. <tim.triche@gmail.com>
}
\references{
  Jaffe AE, Murakami P, Lee H, Leek JT, Fallin MD, Feinberg AP,
  Irizarry RA (2012) Bump hunting to identify differentially
  methylated regions in epigenetic epidemiology studies.
  International Journal of Epidemiology 41(1):200-9.
}
\examples{

## not run:
##
## bumps2bed(bumps.culturedCD34)
## 
## # which will produce a file 'bumps.culturedCD34.bed' by deparsing the call...
## 

}
