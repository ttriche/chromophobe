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

setClass("Segmentation", contains="GRanges", 
        # {{{ a GRanges with defined states
         representation(states="StatesORNULL")
        ) # }}}
Segmentation <- function(from, s=NULL) { # {{{
  new('Segmentation', from, states=s)
} # }}}
setAs("GRanges", "Segmentation", function(from) { # {{{
  Segmentation(from)
}) # }}}

setClass("SegmentationList", contains="GRangesList",
         # {{{ a GRangesList with defined states
         representation(states="StatesORNULL")
        ) # }}}
SegmentationList <- function(from, s=NULL) { # {{{
  if(class(from) == 'list') from <- GRangesList(from)
  new('SegmentationList', endoapply(from, function(x) Segmentation(x)),states=s)
} # }}}
setAs("GRangesList", "SegmentationList", function(from) { # {{{
  SegmentationList(from)
}) # }}}

setClass('JointSegmentation', contains="SummarizedExperiment",
    # {{{ a tweaked SummarizedExperiment
    representation(probinit='dataframeORNULL',
                   emissions='dataframeORNULL', 
                   transitions='matrixORNULL',
                   states='StatesORNULL',
                   rowData='SegmentationList')) # }}}
setAs("JointSegmentation", "SegmentationList", function(from) { # {{{
  SegmentationList(rowData(from), states(from))
}) # }}}
validJointSegmentation <- function(object) { # {{{
  if(identical(rownames(object), colnames(object))) TRUE
  else paste('Non-identical row and column names:',
             paste(setdiff(colnames(object), rownames(object)), collapse=', '),
             'not in colnames, ', 
             paste(setdiff(rownames(object), colnames(object)), collapse=', '),
             'not in rownames!')
} # }}}
# setValidity("JointSegmentation", validJointSegmentation) ## rownames==colnames

setClass('Occupancy',contains="DataFrame",representation(states="StatesORNULL"))
Occupancy <- function(DF, statesData=NULL) { # {{{
  occ <- new('Occupancy')
  for(i in slotNames(DF)) slot(occ, i) <- slot(DF, i)
  slot(occ, 'states') <- statesData
  return(occ)
} # }}}
setAs("DataFrame", "Occupancy", function(from) { # {{{
  Occupancy(from)  
}) # }}}
setAs("Occupancy", "data.frame", function(from) { # {{{
  class(from) <- 'DataFrame'
  as.data.frame(from)
}) # }}}

setClass("TrackHub", ## for downloading HMMs, etc -- may move to rtracklayer
         representation(hubUrl = "character", # {{{
                        hubName = "character",
                        hubNotes = "SimpleList")) # }}}
