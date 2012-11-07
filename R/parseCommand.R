parseCommand <- function(command) { # parse commands from webpage.html
  argv <- list()
  settings <- strsplit(strsplit(command, ':')[[1]][2], ' ')[[1]][-1]
  ## typical: "LearnModel -s 123 -r 25 -x 15000 HMM.binary HMM.model20 20 hg19"
  argv$command <- settings[1]
  if('-x' %in% settings) argv[['maxtime']] <- settings[grep('-x', settings)+1]
  if('-r' %in% settings) argv[['maxiter']] <- settings[grep('-r', settings)+1]
  if('-s' %in% settings) argv[['seed']] <- settings[grep('-s', settings)+1]
  if('-printposterior' %in% settings) argv$posterior <- TRUE
  mandatory <- settings[(length(settings)-3):length(settings)]
  names(mandatory) <- c('binaryDir','modelDir','modelStates','genome')
  for(nm in names(mandatory)) argv[[nm]] <- as.character(mandatory[nm])
  return(argv) # experimentData?
}
