testTrackHub <- FALSE
if(testTrackHub==TRUE) { # {{{
  
require(rtracklayer)
#
# based on the hub track spec (such as it is) from UCSC, found at 
# http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html and 
# http://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html
#
parseStanza <- function(x) { # {{{ parse lines as rowname rowvalue rowvalue...
  if(!is.null(x)) {
    names(x) <- sapply(x, function(y) strsplit(y, ' ')[[1]][1])
    x <- sapply(x, function(y) paste(unlist(strsplit(y,' ')[[1]][-1]),
                                     collapse=' '))
    return(as.list(x))
  } else {
    return(NULL)
  }
} # }}}
#
getStanzas <- function(con) { # {{{ get separated stanzas as a list
  stanzas <- list()
  x <- readLines(con)
  ss <- c()
  for(i in seq_len(length(x))) { 
    x[i] <- gsub('^\ *', '', gsub('\ *$', '', x[i]))
    if(x[i] == '') {
      stanzas[[length(stanzas)+1]] <- parseStanza(ss)
      ss <- c()
    } else {
      ss <- c(ss, x[i])
    }
  }
  # special case: one stanza, no empty lines (not really so special)
  if(length(ss) > 0) stanzas[[length(stanzas)+1]] <- parseStanza(ss)
  return(stanzas)
} # }}}
#
getTrack <- function(tracks) { # {{{ get separated stanzas as a list
  info <- c(track=NA, bigDataUrl=NA, parent=NA)
  for(i in names(info)) if(i %in% names(tracks)) info[i] <- tracks[i]
  return(unlist(info))
} # }}}
#
getTrackDb <- function(trackDbs) { # {{{ get separated stanzas as a list
  db <- as.data.frame(do.call(rbind, lapply(trackDbs, getTrack)))
  rownames(db) <- db$track
  db[, -1]
} # }}}
#
getTrackDbs <- function(stanzas) { # {{{ get separated stanzas as a list
  genomes <- lapply(stanzas, function(x) x[['genome']])
  tracks <- lapply(stanzas, function(x) x[['trackDb']])
  names(tracks) <- genomes
  return(tracks)
} # }}}
#
hubUrl <- 'http://vizhub.wustl.edu/VizHub'
hubFile <- 'RoadmapReleaseAll.txt'
hubSpec <- parseStanza(readLines(url(paste(hubUrl, hubFile, sep='/'))))
##
## hub VizHub
## shortLabel Roadmap Epigenomics Data Complete Collection at Wash U VizHub
## longLabel Roadmap Epigenomics Human Epigenome Atlas Data Complete Collection
## genomesFile genomesall.txt
## email twang@genetics.wustl.edu
##
hubGenomes <- getStanzas(url(paste(hubUrl, hubSpec$genomesFile, sep='/')))
hubGenome <- 'hg19'
hubTrackDbs <- getTrackDbs(hubGenomes)
## 
hubTrackDbUrl <- url(paste(hubUrl, hubTrackDbs[[hubGenome]], sep='/'))
hubTrackDb <- getStanzas(hubTrackDbUrl)
hubTracks <- getTrackDb(hubTrackDb)
hubParents <- unique(hubTracks$parent[which(!is.na(hubTracks$parent))])
names(hubParents) <- hubParents
hubDag <- lapply(hubParents, function(x) hubTracks[which(hubTracks$parent==x),])
hubDag <- lapply(hubDag, function(x) {
  x$dataUrl <- paste(hubUrl, hubGenome, x$bigDataUrl, sep='/')
  return(x)
})
hubDagUrls <- lapply(hubDag, function(x) x[,'dataUrl',drop=F])
head(hubDagUrls)

}  # }}}
