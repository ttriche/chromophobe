require(rtracklayer)

setClass("TrackHubSession",
         representation(url = "character", 
                        hubUrl = "character",
                        hguid = "numeric",
                        views = "environment"),
         contains = "UCSCSession")

setClass("TrackHubQuery",
         representation(session = "TrackHub",
                        track = "characterORNULL",
                        table = "characterORNULL",
                        range = "GRanges",
                        outputType = "characterORNULL",
                        NAMES = "characterORNULL",
                        intersectTrack = "characterORNULL"),
         contains = "UCSCTableQuery")



