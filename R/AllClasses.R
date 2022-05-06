# definition (document me!)
# mostly this is so that S4 methods for plotting etc can be constructed 
# and dispatched appropriately, as well as checkers for necessary bits 
setClass("GenomicSegmentation",
         representation(trackLine = "TrackLine"),
         contains = "GRanges")


#' A subclass of GRanges with a trackLine slot, similar to UCSCData
#' 
#' And also a friendly constructor for such things 
#' 
#' @param   gr    a GRanges or UCSCData object, as from rtracklayer::import.bed
#' @param   ...   additional arguments, currently ignored 
#' 
#' @return        a GenomicSegmentation object, very similar to UCSCData 
#' 
#' @import        GenomicRanges
#' @import        rtracklayer 
#' 
#' @export 
GenomicSegmentation <- function(gr, ...) {
  gs <- as(gr, "GenomicSegmentation")
  if ("trackLine" %in% slotNames(gr)) gs@trackLine <- gr@trackLine
  return(gs)
}


# trackName accessors and setters for common data structures 
setMethod("trackName", "UCSCData", function(x) x@trackLine@name)
setReplaceMethod("trackName", c("UCSCData", "character"),
                 function(x, value) { x@trackLine@name <- value; x })
setMethod("trackName", "GenomicSegmentation", function(x) x@trackLine@name)
setReplaceMethod("trackName", c("GenomicSegmentation", "character"),
                 function(x, value) { x@trackLine@name <- value; x })


# slight tweak of UCSCData `show` method
setMethod("show", "GenomicSegmentation", 
          function (object) {
            if (!is.null(trackName(object))) {
              cat(class(object),", trackName:'",trackName(object),"'\n", sep="")
            }
            callNextMethod()
          })


# subclass for straight non-overlapping HMM
setClass("ChromHMM", contains = "GenomicSegmentation")
setValidity("ChromHMM", function(object) TRUE )


# subclass for stacked HMM derived from >1 HMM tracks
setClass("StackedHMM", contains = "GenomicSegmentation")
setValidity("StackedHMM", function(object) TRUE )


# subclass for scChromHMM data import
setClass("scChromHMM", contains = "GenomicSegmentation")
setValidity("scChromHMM", function(object) TRUE )


# somewhat of a placeholder class for now 
# todo: back and forth to/from RaggedExperiment
# todo: add coercion(s) to/from a segmenter::segmentation object 
setClass("GenomicSegmentationList",
         representation(trackLines = "List"),
         contains = "List")
setValidity("GenomicSegmentationList", function(object) TRUE)


#' Like GRangesList, but with a trackLine slot, similar to UCSCData
#' 
#' And also a friendly constructor for such things 
#' 
#' @param   grl   a GRangesList, with names that become trackLine names
#' @param   ...   additional arguments, currently ignored 
#' 
#' @return        a GenomicSegmentationList object
#' 
#' @import        GenomicRanges
#' @import        rtracklayer 
#' 
#' @export 
GenomicSegmentationList <- function(grl, ...) {
  
  lst <- lapply(grl, as, "GenomicSegmentation")
  for (i in names(grl)) trackName(lst[[i]]) <- i

  # add a validity method to match length(gsl@trackLines) and length(gsl) 
  trackLines <- List(lapply(lst, slot, "trackLine"))
  gsl <- new("GenomicSegmentationList", trackLines=trackLines)
  stop("GenomicSegmentationList coercion is incompletely implemented :-/")

}


# somewhat of a placeholder class for now 
setClass("ChromHMMList", contains = "GenomicSegmentationList")
setValidity("ChromHMMList", function(object) TRUE )


# somewhat of a placeholder class for now 
setClass("StackedHMMList", contains = "GenomicSegmentationList")
setValidity("StackedHMMList", function(object) TRUE )


# scChromHMMList is kind of its own special beast...
