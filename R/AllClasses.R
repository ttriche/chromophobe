setClass('States', contains="DataFrame") 
setClassUnion('StatesOrNULL', c('NULL','States'))
States <- function(DF) { # {{{
  stopifnot(any(c('Id','Name') %in% names(DF)))
  if(!'Name' %in% names(DF)) DF[['Name']] <- DF[['Id']]
  if(!'Id' %in% names(DF)) DF[['Id']] <- DF[['Name']]
  if(!'Color' %in% names(DF)) DF[['Color']] <- colors()[seq_along(DF[['Id']])]
  class(DF) <- 'States'
  return(DF)
} # }}}
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
setAs("States", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}

setClass("Segmentation", contains="GRanges")
setAs("GRanges", "Segmentation", function(from) { # {{{
  class(from) <- 'Segmentation'
  return(from)
}) # }}}

setClass("SegmentationList", contains="GRangesList")
SegmentationList <- function(...) { # {{{
    as(GRangesList(list(...)), 'SegmentationList')
} # }}}
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
    representation(probinit='data.frame',
                   emissions='data.frame', 
                   transitions='matrix',
                   states='StatesOrNULL',
                   rowData='SegmentationList')) # }}}

setClass('Occupancy', contains="DataFrame")
Occupancy <- function(DF) { # {{{
  class(DF) <- 'Occupancy'
  return(DF)
} # }}}
setAs("Occupancy", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}

setClass("TrackHub", 
         representation(hubUrl = "character", # {{{
                        hubName = "character",
                        hubNotes = "SimpleList")) # }}}

setClass("TrackHubQuery", contains = "UCSCTableQuery")

setClass("TrackHubSession", contains = c("UCSCSession","TrackHub")) 
TrackHubSession <- function(hubUrl, ... ) { # {{{
  new('TrackHubSession', hubUrl=hubUrl, ...)
} # }}}
