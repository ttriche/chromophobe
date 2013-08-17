getSNPplots<-function(x, individuals=NULL, rotate=FALSE, col.fun='SNP') {
  require('GMD')
  require('fastcluster')
  tmp <- matrix()
  if(is(x, 'SummarizedExperiment')) {
    probes <- grep('^rs', names(rowData(x)))
    snps <- grep('^rs', names(rowData(x)), val=T)
    tmp <- assays(x, withDimnames=F)[[1]][ probes, ]
  } 
  if(nrow(tmp) < 2) {
    stop('Need a SummarizedExperiment with SNP (rsXX) features... none found')
  }
  tmp = round(tmp*2)
  if(is.null(individuals)) individuals = dim(x)[2]
  SNP <- colorRampPalette(c('blue','yellow','red'))
  bw <- colorRampPalette(c('white','gray','black'))
  heading <- paste('SNPs for', individuals, 'individuals')
  if(rotate) {
    capture.output({ # {{{
      clusts <- suppressWarnings(heatmap.3(t(tmp), scale="none", trace="none", 
                                 color.FUN=get(col.fun), dendrogram='none', 
                                 labCol=snps, kr=individuals, Colv=T, Rowv=T, 
                                 labRow=colnames(x), main=heading))
      clusts$clusters <- clusts$row.clusters
      clusts$ind <- clusts$rowInd
    }) # }}}
  } else { 
    capture.output({ # {{{
      clusts <- suppressWarnings(heatmap.3(tmp, scale="none", trace="none", 
                                 color.FUN=get(col.fun), dendrogram='none', 
                                 labRow=snps, kc=individuals, Colv=T, Rowv=T, 
                                 labCol=colnames(x), main=heading))
      clusts$clusters <- clusts$col.clusters
      clusts$ind <- clusts$colInd
    }) # }}}
  }
  individual <- clusts$clusters
  individual[clusts$ind] <- individual
  individual <- as.factor(individual)
  message('Assigned identity for each sample:')
  return(individual)
}
