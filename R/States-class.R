setClass('States', contains="DataFrame") ## yep, another one of these...
States <- function(DF) { # {{{
  stopifnot(any(c('Id','Name') %in% names(DF)))
  if(!'Name' %in% names(DF)) DF[['Name']] <- DF[['Id']]
  if(!'Id' %in% names(DF)) DF[['Id']] <- DF[['Name']]
  if(!'Color' %in% names(DF)) DF[['Color']] <- colors()[seq_along(DF[['Id']])]
  class(DF) <- 'States'
  return(DF)
} # }}}
setAs("States", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}
setMethod('$', 'States', # {{{
  function(x, name) {
    xx <- x[[name, exact=FALSE]]
    names(xx) <- rownames(x)
    return(xx)
  }) # }}}

## color name to RGB, NA if it fails
toRgb <- function(color) { # {{{
  if(length(color) > 1) { 
    res <- t(sapply(color, toRgb))
    colnames(res) <- c('red','green','blue')
    return(res)
  } else {
    x <- try(col2rgb(color))
    if(!inherits(x, "try-error")) return(x)
    else return(NA)
  }
} # }}}

## color name or RGB vector to hex, NA if it fails 
toHex <- function(RGB) { # {{{
  if(is(RGB,'matrix')) {
    res <- apply(RGB, 1, toHex)
    names(res) <- rownames(RGB)
    return(res)
  } else if(is(RGB, 'character')) {
    return(toHex(toRgb(RGB)))
  } else if(length(RGB) == 3 && is.integer(RGB)) {
    return(rgb(red=RGB[1], green=RGB[2], blue=RGB[3], max=255))
  } else {
    return(NA)
  }
} # }}}

## is everything set up so that it won't blow to pieces?
setValidity("States", function(object) { # {{{
  msg <- NULL
  for(i in c('Id','Name','Color')) if(!i %in% names(object)) {
     msg <- validMsg(msg, sprintf("object of class '%s' needs column '%s'", 
                                  class(object), i))
  }
  if(!identical(rownames(object), object[['Id']])) {
     msg <- validMsg(msg, "rownames(object) do not match object$Id")
  }
  if(any(is.na(toRgb(object[['Color']]))) && !all(is.na(object[['Color']]))) {
     msg <- validMsg(msg, "object$Color contains invalid, non-NA R color names")
  }
  if (is.null(msg)) TRUE else msg
}) # }}}

## various accessors to make my life easier
setGeneric('states', function(object, ...) standardGeneric('states'))
setMethod('states', signature(object='States'), # {{{
          function(object) states(object$Names)) # }}}

setGeneric('stateIds', function(object, ...) standardGeneric('stateIds'))
setMethod('stateIds', signature(object='States'), # {{{
          function(object) object$Id) # }}}

setGeneric('stateNames', function(object, ...) standardGeneric('stateNames'))
setMethod('stateNames',  signature(object='States'), # {{{
          function(object) object$Name) # }}}

setGeneric('stateColors', function(object, ...) standardGeneric('stateColors'))
setMethod('stateColors', signature(object='States'), # {{{
          function(object)object$Color) # }}}

setGeneric('stateRGB', function(object, ...) standardGeneric('stateRGB'))
setMethod('stateRGB', signature(object='States'), # {{{
          function(object) toRgb(stateColors(object))) # }}}

setGeneric('stateHex', function(object, ...) standardGeneric('stateHex'))
setMethod('stateHex', signature(object='States'), # {{{
          function(object) toHex(stateRGB(object))) # }}}
