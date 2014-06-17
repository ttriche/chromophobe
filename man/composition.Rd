\name{composition}
\alias{composition}
\title{
  Compositional tools: estimating cell counts & using simplex transformations
}
\description{
  John Aitchison has repeatedly argued for additive log-contrasts as the proper
  form in which to manipulate compositional data (that is, data where the parts
  must sum to unity, and thus induce negative correlation as an artifact).  It 
  turns out that simplicial principal components are also a terrific way to 
  capture cell-type-related artifactual effects and control for them in DNA 
  methylation and chromatin immunopreciptation studies.  The utility functions
  in chromophobe make this a straightforward practice to implement.
}
\details{
  See the 'compositions' package and accompanying book for much, much more.
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
