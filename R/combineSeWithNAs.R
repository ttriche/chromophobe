## Take a bunch of ragged SummarizedExperiment objects and inject NAs 
## for the ones which don't have observations at a feature where others do.
## This reduces to "find the superset of features and inject NAs elsewhere".
## The assumption is that some sort of smoothing will be done later on. 
## This was originally written to combine enhanced RRBS data for comparisons.
##
combineSeWithNAs <- function(x, y, ...) { 

  eq <- findOverlaps(rowData(x), rowData(y), type='equal')
  ol <- findOverlaps(rowData(x), rowData(y), type='any')
  if(length(ol) != length(eq)) {
    message('Some features overlap only partially -- you may want to disjoin()')
  }

  ## would prefer to automatically do this if x is itself a list!
  if(length(list(...)) > 0L) y = do.call(combineSeWithNAs, list(y, ...))

  ## get some annoying preliminaries out of the way first...
  if(class(x)!=class(y)) stop(paste("class mismatch: x is a", class(x), ", ", 
                                     "while y is a", class(y), sep=""))
  if(names(colData(x)) != names(colData(y))) 
    stop("names(colData()) differ between x and y")
  if((ncol(values(x)) > 0 || ncol(values(y)) > 0) &&
     (names(values(x)) != names(values(y))))
    stop("names(values(rowData())) differ between x and y")
  if(names(assays(x, withDimnames=F)) != names(assays(y, withDimnames=F))) 
    stop("names(assays()) differ between x and y")

  ## superset of features/row data in common
  asy.names = names(assays(x, withDimnames=F))
  row.dat = sort(union(rowData(x), rowData(y)))
  message("Adding NAs for unshared features...")
  asy.dat = lapply(assays(x, withDimnames=F), function(x) 
                     matrix(NA, ncol=ncol(x)+ncol(y), nrow=length(row.dat)))

  ## ragged merge of x
  xcols = seq_len(ncol(x))
  o = findOverlaps(row.dat, rowData(x))
  for(i in asy.names) 
    asy.dat[[i]][queryHits(o),xcols] = assays(x,withDim=F)[[i]][subjectHits(o),]

  ## ragged merge of y
  ycols = ncol(x) + seq_len(ncol(y))
  o = findOverlaps(row.dat, rowData(y))
  for(i in asy.names)
    asy.dat[[i]][queryHits(o),ycols] = assays(y,withDim=F)[[i]][subjectHits(o),]

  ## merge column data if possible
  col.dat <- merge(colData(x), colData(y), all=TRUE)

  ## ugh
  gc(,T) 

  ## return the merged SummarizedExperiment
  SummarizedExperiment(assays=asy.dat, rowData=row.dat, colData=col.dat)

}
