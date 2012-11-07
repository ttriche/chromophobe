importBinarized <- function(HMM, binarizedPath) {

  if(!is(HMM, 'JointSegmentation')) {
    stop('Need an already-loaded JointSegmentation model to load binary files!')
  }

  ## given an HMM GRanges/GRangesList, read in the matrices of marks (0/1/NA)
  ## corresponding to the binned/binarized data from ChromHMM BinarizeBed and
  ## stuff them into sparse Matrix objects in a SummarizedExperiment
  ##
  ## This may require some serious trickery in order to be doable properly in R
  ## 
  stop("Ain't done yet.  Feel free to send me a patch though! --t")

}

