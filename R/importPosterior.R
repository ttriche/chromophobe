import.posterior <- function(HMM, posteriorPath, states=NULL) {

  ## given an HMM GRanges/GRangesList, turn it into a SummarizedExperiment, by
  ## loading the posterior probabilities into a sparse Matrix for each GRanges
  ## (so if the HMMs are in a GRangesList, the result will be a list of SEs)
  ##
  ## FIXME: make sure to distinguish a list of GRanges from a GRangesList, 
  ##        in the sense that a GRangesList with full genome coverage is a very
  ##        different animal from a per-state GRangesList that can be requested
  ##        from import.ChromHMM for producing pyro plots.  The latter must be 
  ##        unlisted back into a single GRanges in order to cough up an HMM SE.
  ## 
  stop("Ain't done yet.  Feel free to send me a patch though! --t")

}

