## for JointSegmentation assembly, e.g. Roadmap
setClassUnion('matrixORNULL', c('NULL','matrix'))
setClassUnion('dataframeORNULL', c('NULL','data.frame'))

setClass('States', contains="DataFrame") 
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
setAs("DataFrame", "States", function(from) { # {{{
  States(from)  
}) # }}}
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
setClassUnion('StatesORNULL', c('NULL','States'))

setClass("Segmentation", contains="GRanges")
setAs("GRanges", "Segmentation", function(from) { # {{{
  class(from) <- 'Segmentation'
  return(from)
}) # }}}

setClass("SegmentationList", contains="GRangesList")
setAs("GRangesList", "SegmentationList", function(from) { # {{{
  from <- endoapply(from, function(x) {
    class(x) <- 'Segmentation'
    return(x)    
  }) 
  class(from) <- 'SegmentationList'
  return(from)
}) # }}}

setClass('JointSegmentation', contains="SummarizedExperiment",
    # {{{ a tweaked SummarizedExperiment
    representation(probinit='dataframeORNULL',
                   emissions='dataframeORNULL', 
                   transitions='matrixORNULL',
                   states='StatesORNULL',
                   rowData='SegmentationList')) # }}}

setClass('Occupancy',contains="DataFrame",representation(states="StatesORNULL"))
Occupancy <- function(DF, statesData=NULL) { # {{{
  class(DF) <- 'Occupancy'
  DF@states <- statesData
  return(DF)
} # }}}
setAs("DataFrame", "Occupancy", function(from) { # {{{
  Occupancy(from)  
}) # }}}
setAs("Occupancy", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}

setClass("TrackHub", 
         representation(hubUrl = "character", # {{{
                        hubName = "character",
                        hubNotes = "SimpleList")) # }}}
