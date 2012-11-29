require(rtracklayer)

setClass("TrackHub", 
         representation(hubUrl = "character",
                        hubName = "character",
                        hubNotes = "SimpleList"))

setClass("TrackHubSession",
         contains = c("UCSCSession","TrackHub"))

setClass("TrackHubQuery",
         representation(session = "TrackHub",
                        range = "GRanges",
                        outputType = "characterORNULL",
                        NAMES = "characterORNULL",
                        intersectTrack = "characterORNULL"),
         contains = "UCSCTableQuery")



