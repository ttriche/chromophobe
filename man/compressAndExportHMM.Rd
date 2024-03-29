% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compressAndExportHMM.R
\name{compressAndExportHMM}
\alias{compressAndExportHMM}
\title{Reduce a (usually simplified) segmentation, export to BED, compress, index.}
\usage{
compressAndExportHMM(HMM, name, filename = NULL)
}
\arguments{
\item{HMM}{a GenomicSegmentation or GRanges object}

\item{name}{a name to include as a trackLine in the BED file}

\item{filename}{the output (BED) file path (default: ./name.genome.bed)}
}
\description{
Reduce a (usually simplified) segmentation, export to BED, compress, index.
}
\examples{

data(chr19_HMM, package="chromophobe")
data(remc18state, package="chromophobe")

simpler <- remc18state
simpler$SIMPLE <- sub("(Promoter|Enhancer)", "Active", simpler$SIMPLE)
simpler[simpler$SIMPLE == "Active", "RGBSIMPLE"] <- "255,0,0" # Red
simpler$SIMPLE <- sub("(Transcribed|Het_Rpt_Qui)", "Other", simpler$SIMPLE)
simpler[simpler$SIMPLE == "Other", "RGBSIMPLE"] <- "255,255,255" # White

tmpdir <- tempdir()
setwd(tmpdir)
simplerHMM <- simplifyHMM(chr19_HMM, cols=simpler)
compressAndExportHMM(simplerHMM, name="chr19_HMM_simple", 
                     filename=file.path(tmpdir, "chr19_HMM.simple.mm10.bed"))
list.files(tmpdir) 

}
