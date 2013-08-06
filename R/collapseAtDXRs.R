## Collapse a SummarizedExperiment-like object over a GRanges of DXRs, sensibly.
## Handy for things like chromophobe model testing, GenometriCorr, lars, etc.
## 
## arguments:
##  x           ==  anything descended from a SummarizedExperiment
##  y           ==  differentially whatevered regions (DXRs) as a GRanges
##  how         ==  how to summarize the assay values (default: median)
##  parallel    ==  run the assays in parallel?  (default: false)
##
## value:
##  xy          ==  an object with same colData but with new rowData and assays
## 
## FIXME: let the user provide their own summarizer function
## 
collapseAtDXRs <- function(x, y, how=c('median','mean','sum','max','min'), parallel=F) {  # {{{ 

    ## check arguments
    require(matrixStats)
    require(GenomicRanges)
    stopifnot(is(y, 'GenomicRanges'))
    stopifnot(is(x, 'SummarizedExperiment'))
    xx <- subsetByOverlaps(x, y)

    ## find and index the runs
    hitz <- findOverlaps(xx, y)
    byDXR <- seqapply(split(hitz, subjectHits(hitz)), queryHits)
    names(byDXR) <- paste0('region', length(byDXR))

    ## obtain the summarizer
    how <- match.arg(how)
    fnBy <- switch(how,
                   'median'='colMedians',
                   'mean'='colMeans',
                   'sum'='colSums',
                   'max'='colMaxs',
                   'min'='colMins')

    ## loop through the assays as found in the SE 
    summarizeAsy <- function(xxy) collapseMatrix(xxy, byDXR, fnBy)
    summarizedAssays <- asyApply(xx, summarizeAsy, parallel=parallel)
    
    ## construct a clone of the original SE but with the new correct # of rows 
    res <- x[ seq_along(y), ] ## kludge to retain everything while re-sizing 
    for(i in names(summarizedAssays)) assays(res)[[i]] <- summarizedAssays[[i]]
    rownames(res) <- paste0(names(byDXR), '.', how)

    ## now fix the rowData for the result
    ranges(rowData(res)) <- ranges(y)
    seqdummy <- as(seqnames(res), 'factor') ## Herve fix
    seqnames(rowData(res)) <- seqdummy[as.numeric(match(seqnames(y), seqdummy))]
    strand(rowData(res)) <- strand(y)
    return(res)

} # }}}


## collapse variable, and perhaps overlapping, rows of a matrix-like structure 
##
## arguments:
##  x           ==  anything descended from a SummarizedExperiment
##  z           ==  the listed matrix indices within DXRs (NOT THE ACTUAL DXRs!)
##  fn          ==  the name of the function to be used to summarize x over z 
##
## value:
##  asys        ==  a list of assays, having been transformed
## 
collapseMatrix <- function(x, z, fn) { # {{{
    fn <- selectMethod(fn, class(x))
    do.call(rbind, lapply(z, function(w) fn(x[w,,drop=F])))
} # }}}


## apply a function to each of the assays in an SE, efficiently
##
## arguments:
##  x           ==  anything descended from a SummarizedExperiment
##  fn          ==  the actual function to run on each assay matrix
##  parallel    ==  run the assays in parallel?  (default: false)
##
## value:
##  asys        ==  a list of assays having been transformed
## 
asyApply <- function(x, fn, ..., parallel=FALSE) { # {{{ 
    if(parallel) {
        mclapply(assays(x, withDimnames=FALSE), fn, ...)
    } else { 
        lapply(assays(x, withDimnames=FALSE), fn, ...)
    }
} # }}}
