getXplots<-function(x, nFeats=500, col.fun='jet', pval=.000001) {

  ## FIXME: get rid of dependence on GMD
  require('GMD')
  require('matrixStats')
  require('fastcluster')
  require('GenomicRanges')
  if(!is(x, 'SummarizedExperiment')) {
    stop('You need to provide a SummarizedExperiment for this to work')
  } else if( !('chrX' %in% runValue(seqnames(x))) ){
    stop('Your SummarizedExperiment does not contain any chrX probes!')
  }

  ## order by SD, then select top nFeats CpG loci 
  x.chrX.SD <- rowSds(assays(keepSeqlevels(x, 'chrX'))[[1]])
  Xloci <- names(head(sort(x.chrX.SD, decreasing=TRUE), nFeats))
  XInd <- which(rownames(x) %in% grep('^cg', Xloci, val=T))

  ## FIXME: when CN data is available (e.g. hm450 chips), that may be preferable
  if(!is.null(colData(x)$gender)) {
    colSide <- ifelse(is.na(colData(x)$gender), 'white',
                      ifelse(substr(toupper(colData(x)$gender),1,1) == 'M', 
                             'lightblue', 'pink'))
  } else {
    colSide = rep('white', ncol(x))
  }
  name <- as.character(match.call()["x"])
  asy <- names(assays(x, withDimnames=F))[[1]]
  tmp <- assays(x, withDimnames=F)[[asy]][XInd, ]
  jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                            "yellow", "#FF7F00", "red", "#7F0000"))
  bw <- colorRampPalette(c('white','gray','black'))
  capture.output({
        clusts <- suppressWarnings(
                    heatmap.3(tmp, scale="none", trace="none", 
                              ColIndividualColors=colSide, 
                              color.FUN=get(col.fun), 
                              dendrogram='none', kc=2, labCol=colnames(x), 
                              labRow=Xloci, Colv=T, Rowv=TRUE,
                              main=paste('chrX clustering for', name)))
  })

  by.clust <- lapply(unique(clusts$col.clusters), 
                     function(x) {
                       cm <- colMeans(tmp[, clusts$colInd], na.rm=TRUE)
                       cm[clusts$col.clusters == x]
                     })

  ## do a test to quantify the chance they're all the same
  p <- wilcox.test(by.clust[[1]], by.clust[[2]])$p.value
  if(p > pval) {
    message('Cannot distinguish sex: Pr(one sex) = ',p,' > ',pval,' (cutoff)')
    capture.output({
      clusts <- suppressWarnings(
                  heatmap.3(tmp, scale="none", trace="none", 
                            ColIndividualColors=colSide, 
                            color.FUN=get(col.fun),
                            dendrogram='none', labCol=colnames(x), 
                            labRow=Xloci, Colv=T, Rowv=TRUE,
                            main=paste('chrX clustering for', name)))
    })
    return(rep(NA, ncol(tmp)))
  } else { 
    message('Pr(all samples are the same sex) = ', p)
    labels <- c('M','F')
    ## I feel dirty about this
    sex <- clusts$col.clusters
    sex[clusts$colInd] <- sex
    ## higher == female (XX with one Xi) 
    means <- do.call(c, lapply(by.clust, mean, na.rm=TRUE))
    if(means[1] > means[2]) labels <- rev(labels)
    message('Assigned chrX cluster (assuming this is DNA methylation data):')
    return(gsub('1', labels[1], gsub('2', labels[2], sex)))
  }
}
