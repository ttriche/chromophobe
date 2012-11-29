getVignetteFromSvn <- function(pkg,subDir='vignettes',vignette=NULL,run=FALSE) {

  if(is.null(vignette)) vignette <- pkg # could be different

  ## bookkeeping
  userPwd <- 'readonly:readonly'
  vignetteFile <- paste(vignette, 'Rnw', sep='.')
  stangleFile <- sub('Rnw','R',vignetteFile)
  svnUrl <- 'https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks'
  completeUrl <- paste(svnUrl, pkg, subDir, vignetteFile, sep='/')

  ## the real work 
  require(RCurl)
  Rnw <- getURL(completeUrl, userpwd=userPwd)
  cat(Rnw, file=vignetteFile)
  Stangle(vignetteFile)
  ## Writing to file easyRNASeq.R 

  ## now run it, if you wish 
  if(run==TRUE) source(stangleFile, verbose=T)

}
