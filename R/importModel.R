importModel <- function(file, states=NULL, loud=FALSE) {

  if(loud) message(paste('Importing model parameters from', file, '...'))

  require(utils)
  require(reshape2)
  w <- rep('character', 6)
  model <- suppressWarnings(read.table(file, header=F, skip=1, fill=T, colCl=w))
  names(model)[1] <- 'what'
  model$what <- as.factor(model$what)
  model.pieces <- split(model, model$what)

  ## handy for when we don't have state names 
  getStates <- function(x) unique(as.numeric(x[['state']]))

  ## initial state probabilities
  probinit <- model.pieces$probinit[,2:3]
  names(probinit) <- c('state','p')
  if(is.null(states)) {
    Id <- getStates(probinit)
    Name <- paste0('E', Id)
    states <- States(DataFrame(Id, Name))
  } else {
    message('Collapsing model states has not been debugged -- you are warned!')
  }
  probinit$state <- with(probinit,factor(stateNames(states)[as.numeric(state)]))
  probinit$p <- as.numeric(probinit$p) 
  model.pieces$probinit <- probinit

  ## emission probabilities
  emissions <- model.pieces$emissionprobs[,-1]
  names(emissions) <- c('state','id','mark','keep','p')
  emissions$state <- with(emissions, 
                          factor(stateNames(states)[as.numeric(state)]))
  emissions$mark <- with(emissions, factor(mark, labels=unique(mark)))
  emissions$p <- as.numeric(emissions$p)
  emissions <- emissions[ emissions$keep==1, c('state','mark','p') ]
  model.pieces$emissions <- emissions

  ## transition probabilities
  transitions <- model.pieces$transitionprobs[,-1]
  names(transitions) <- c('from','to','p')
  transitions$from <- as.numeric(transitions$from)
  transitions$to <- as.numeric(transitions$to)
  transitions$p <- as.numeric(transitions$p)
  transitions <- transitions[, c('from','to','p')]
  model.pieces$transitions <- acast(transitions, to ~ from, value.var='p')

  ## model states 
  model.pieces$states <- states

  return(model.pieces[ c('probinit','emissions','transitions','states')])

}
