import.model <- function(file, states=NULL) {

  require(utils)
  w <- rep('character', 6)
  model <- suppressWarnings(read.table(file, header=F, skip=1, fill=T, colCl=w))
  names(model)[1] <- 'what'
  model$what <- as.factor(model$what)
  model.bits <- split(model, model$what)

  emissions <- model.bits$emissionprobs[,-1]
  names(emissions) <- c('state','id','mark','keep','p')
  emissions$state <- as.factor(paste0('E', emissions$state))
  emissions$mark <- as.factor(emissions$mark)
  emissions$p <- as.numeric(emissions$p)
  emissions <- emissions[ emissions$keep==1, c('state','mark','p') ]
  rownames(emissions) <- paste(emissions$state, emissions$mark, sep='.')

  transitions <- model.bits$transitionprobs[,-1]
  names(transitions) <- c('from','to','p')
  transitions$from <- as.factor(paste0('E', transitions$from))
  transitions$to <- as.factor(paste0('E', transitions$to))
  transitions$p <- as.numeric(transitions$p)
  transitions <- transitions[, c('from','to','p')]
  rownames(transitions) <- paste(transitions$from, transitions$to, sep='.')

  model <- list(emissions=emissions, transitions=transitions)
  return(model)

}
