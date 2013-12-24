\name{composition}
\alias{composition}
\docType{function}
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
  Compositional data analysis in R.  (authors)
}
\examples{

## not run:
##
## require(FlowSorted.Blood.450k)
## counts <- getBloodCellCounts(your.GenomicRatioSet.from.GEO)
## plotCellCounts(counts)
## plotCellCountPCs(counts)
## countPCs <- getCellCountPCs(counts, k=2)
## 

}
