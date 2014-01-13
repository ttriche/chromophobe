## Convenience functions for doing "proper" compositional correction of samples
## document in chromophobe manpages & package description!
##

## adapt minfi's cell count estimation to work without IDATs
## also, fix a bug where it can return counts < 0, and normalize
## (Houseman's original code was less convoluted than that in minfi)
getBloodCellCounts <- function(grSet, referenceMset=NULL) { # {{{
  require(minfi)
  require(quadprog)
  require(FlowSorted.Blood.450k)
  if(is.null(referenceMset)) { # {{{
    cat("[estimateCellCounts] Loading reference data FlowSorted.Blood.450k...")
    referenceMset <- get('FlowSorted.Blood.450k')
    referenceMset <- preprocessQuantile(referenceMset)
  } # }}}
  cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran")
  compData <- minfi:::pickCompProbes(referenceMset, cellTypes=cellTypes)
  coefs <- compData$coefEsts
  coefs <- coefs[ intersect(rownames(grSet), rownames(coefs)), ]
  rm(referenceMset)
  cat("[estimateCellCounts] Estimating composition.\n")
  counts <- minfi:::projectCellType(getBeta(grSet)[rownames(coefs), ], coefs)
  rownames(counts) <- sampleNames(grSet)
  if(any(counts < 0)) { ## this happens
    counts[ which(counts < 0) ] <- 0
    counts <- apply(counts, 1, function(x) x/sum(x))
  }
  return(counts) ## how hard was that, really?! 
} # }}}

## stacked bars (for raw counts; note 'compositions' features additional plots)
plotCellCounts <- function(estimates) { # {{{
  require(reshape2)
  if(is(estimates, 'matrix')) {
    estimates <- melt(estimates)
    names(estimates) <- c('subject','celltype','fraction')
  } else { 
    stop("Need matrix of counts with columns 'subject', 'celltype', 'fraction'")
  }
  require(ggplot2)
  p <- ggplot(estimates, aes(y=fraction, x=factor(subject), fill=celltype)) +
              geom_bar(position='fill') + xlab('subject') + ylab('fraction')
  p <- p + scale_fill_brewer(type="div", palette=7) +
           ggtitle('Estimated leukocyte composition of each sample') +
           theme(panel.background=element_blank(),
                 axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=9))
  p
} # }}}

## get simplicial PCs for counts (useful for correcting effects of diff. counts)
getCellCountPCs <- function(estimates, k=2) { # {{{
  require(compositions)
  pr <- princomp(acomp(estimates))
  pr$scores[,1:k]
} # }}}

## visualize the fit and characteristics of closed simplex PCs for cell counts
plotCellCountPCs <- function(estimates, type="scale") { #{{{
  require(compositions)
  pr <- princomp(acomp(estimates))
  dots <- rep(".", times=nrow(estimates))
  plot(pr, xlabs=dots, type="biplot", scale=ifelse(type=='scale', 1, 0))
} # }}}

## e.g.
##
## if(FALSE) {
##   library(minfi)
##   HuGeF.pdat <- readRDS("HuGeF.pData.rds")
##   HuGeF <- read.450k.exp(base=".", targets=HuGeF.pdat)
##   for( i in colnames(HuGeF.counts) ) pData(HuGeF)[,i] <- HuGeF.counts[,i]
##   HuGeF.counts <- estimateCellCounts(HuGeF)
##   HuGeF <- preprocessQuantile(HuGeF)
##   saveRDS(HuGeF, file="HuGeF.rds")
##   HuGeF.PCs <- compPcs(HuGeF.counts)
##   for(i in colnames(HuGeF.PCs)) colData(HuGeF)[,i] <- HuGeF.PCs[,i]
##   saveRDS(HuGeF, file='HuGeF.rds')
##   dmat <- with(as(colData(HuGeF), 'data.frame'),
##                model.matrix(~ age + gender + Comp.1 + Comp.2))
##   stopifnot(nrow(dmat) == ncol(HuGeF)) ## double check that it will work
## 
##   ## HuGeF age DMRs
##   library(doMC)
##   registerDoMC(2) ## for parallel bump hunting ... be careful though
##   HuGeF.age.DMRs <- bumphunter(HuGeF, dmat, pickCutoff=T, cutoffQ=0.99)
##   saveRDS(HuGeF.age.DMRs, file='HuGeF.age.DMRs.rds')
## }

