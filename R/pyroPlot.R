## SE is a SummarizedExperiment with GRanges as rowData, covariates as colData
## GR[L] is a partitioning of the genome which may or may not match rowData
## groupBy is an optional logical variable by which to stratify the tumors 
## assay is if you aren't using the first assay in the SE for some reason
## withNormals specifies whether to look for normals in the SE to split off
## tidy indicates whether SNPs, etc. should be removed before plotting
## 
pyroPlot <- function(SE, GRL, groupBy=NULL, assay=NULL, withNormals=T, tidy=F) {

  require(regulatoR)

  ## toss out bogus probes if required
  ## dropUnclusterables() is from regulatoR
  if(tidy==TRUE) { # {{{
    require(FDb.FDb.UCSC.snp135common.hg19)
    SNPs <- features(FDb.FDb.UCSC.snp135common.hg19)
    SE <- dropUnclusterables(SE, SNPs=commonSNPs) 
  } # }}}

  if(withNormals==TRUE) { # {{{
    SEs <- list(Normal=SE[ , which(SE$histology == 'Normal Tissue')],
                Tumor=SE[ ,  which(SE$histology != 'Normal Tissue')]) # }}}
  } else { # {{{
    SEs = list(Tumor=SE)
  } # }}}
  if(!is.null(groupBy)) { # {{{
    SEs <- lapply(SEs, function(x) 
                  list(Mutant=x[, which(x[[groupBy]]==TRUE)],
                       Wildtype=x[, -which(x[[groupBy]]==TRUE)])) # }}}
  } else { # {{{
    SEs <- lapply(SEs, function(x) list(wildtype=x))
  } # }}}

  ## we end up with SEs = list(Normal=list('xxx mutant'=SE, 'xxx wildtype'=SE),
  ##                            Tumor=list('xxx mutant'=SE, 'xxx wildtype'=SE) )
  ## or just {Tumor, Normal}{wildtype, mutant} with empties. 

  ## reorder samples in each group based on the mean value across features
  browser()
  SEs <- lapply(SEs, function(x) 
                lapply(x, function(y) y[,order(colMeans(asy.fast(y, assay)))]))

  ## should I be permitting multiple GRLs as a list to loop over?
  browser()
  GRL <- subsetByOverlaps(GRL, rowData(SE))
  splitSEs <- lapply(SEs, function(x) lapply(x, function(y) byList(y, GRL)))






  ## now pull out the betas and plot them, split, by status: N, T-wt, T-mut

  # normal
  betaN <- BLCA_beta[,NormID]
  ordN <- names(sort(apply(betaN,2, mean, na.rm=TRUE)))
  # wldtyoe
  betaW <- BLCA_beta[,WtID]
  ordW <- names(sort(apply(betaW,2, mean,na.rm=TRUE)))
  # mutant
  betaM <- BLCA_beta[,MutID]
  ordM <- names(sort(apply(betaM,2, mean,na.rm=TRUE)))
  reorder <- c(ordN, ordW, ordM)
  BLCA_beta <- BLCA_beta[,reorder]

  # pu - promoter upstream 1,500
  # pd  - promoter downstream 1,500
  # p - promoter upstream and dowmstream combined
  # ex - exon
  # ir - intron
  # inter - intergenic

  table(probe.annot$summaryA,exclude=NULL)
  # 
  #         cgi_p        cgi_ex        cgi_ir     cgi_inter     non.cgi_p 
  #         96992         10849         12884         18623         97804 
  #    non.cgi_ex    non.cgi_ir non.cgi_inter  multipleHits          <NA> 
  #         29151         92882         82462          3163             0 
  table(probe.annot$summaryB,exclude=NULL)
  # 
  #        cgi_pu        cgi_pd        cgi_ex        cgi_ir     cgi_inter 
  #         44690         36776         10849         12884         18623 
  #    non.cgi_pu    non.cgi_pd    non.cgi_ex    non.cgi_ir non.cgi_inter 
  #         61231         26889         29151         92882         82462 
  #  multipleHits          <NA> 
  #         28373             0 


  # use summaryA in which promoter upstream and downstream are combined
  # define ESC-H3K27 probes 
  # remove probes that are mapped in multiple categories
  probe.annot <- subset(probe.annot,summaryA!="multipleHits",select=c("summaryA","wgEncodeBroadHistoneH1hescH3k27me3StdPk_0.01"))
  probe.annot$summaryA <- as.character(probe.annot$summaryA)
  cgi_p.K27 <- with(probe.annot, which(summaryA =="cgi_p" & wgEncodeBroadHistoneH1hescH3k27me3StdPk_0.01==1))
  probe.annot$summaryA[cgi_p.K27] <- "cgi_p.PRC" 
  probe.annot <- subset(probe.annot,select=1,drop=FALSE)
  probe.annot$summaryA <- factor(probe.annot$summaryA, 
      levels=c("cgi_p","cgi_p.PRC","non.cgi_p",
          "cgi_ex","non.cgi_ex",
          "cgi_ir", "non.cgi_ir",
          "cgi_inter","non.cgi_inter"))
  table(probe.annot$summaryA)
  # 
  #         cgi_p     cgi_p.PRC     non.cgi_p        cgi_ex    non.cgi_ex 
  #         64100         32892         97804         10849         29151 
  #        cgi_ir    non.cgi_ir     cgi_inter non.cgi_inter 
  #         12884         92882         18623         82462 



  # combine beta values and probe annotation 
  BLCA_beta <- BLCA_beta[rownames(probe.annot),]
  identical(rownames(BLCA_beta),rownames(probe.annot))
  # [1] TRUE

  dat <- cbind(probe.annot, BLCA_beta)

  # remove NA data points
  dat <- na.omit(dat)
  dim(dat)
  # [1] 421785     66


  # split the dataset into different categories  
  dat.split <- split(dat, dat$summaryA)
  dat.split <- lapply(dat.split, function(x) x[,-1])

  # sort probes based on beta values for each sample in each probe category
  meth_heatmap_1 <- function(dat.split){
    dat <- vector("list",length(dat.split))
    names(dat) <- names(dat.split)
    for (i in 1:length(dat.split)){
      temp.dat <- as.matrix(dat.split[[i]])
      for (j in 1:ncol(temp.dat)){
        temp.dat[,j] <- as.numeric(rev(sort(temp.dat[,j])))}
      dat[[i]] <- temp.dat
      rm(temp.dat)}
    dat.split.sorted <<- dat}


  meth_heatmap_1(dat.split) 
  dim.dat <- lapply(dat.split.sorted,dim)
  str(dim.dat)
  # List of 9
  #  $ cgi_p        : int [1:2] 61917 65
  #  $ cgi_p.PRC    : int [1:2] 31880 65
  #  $ non.cgi_p    : int [1:2] 93615 65
  #  $ cgi_ex       : int [1:2] 10383 65
  #  $ non.cgi_ex   : int [1:2] 27515 65
  #  $ cgi_ir       : int [1:2] 12187 65
  #  $ non.cgi_ir   : int [1:2] 88137 65
  #  $ cgi_inter    : int [1:2] 17819 65
  #  $ non.cgi_inter: int [1:2] 78332 65



  ## generate a heatmap
  # Row labels
  library(stringr)
  textlab <- vector(mode="character", length=length(dat.split.sorted))
  for (i in 1:length(textlab))
    textlab[i] <- paste(names(dim.dat)[i],"\n(", str_sub(dim.dat[[i]][1],1,2),",",str_sub(dim.dat[[i]][1],3,5),")",sep="")

  library("matlab")
  pal <- jet.colors(256) 
  breaks <- seq(0,1,length.out=257) 

  meth_heatmap_2 <- function(ls_dat,bot=0.22,lef=5.1,top=0.22,rig=3.1,plotid,text1){
    par(mar = c(bot,lef,top,rig))
    image(x=1:ncol(ls_dat[[plotid]]),y=1:nrow(ls_dat[[plotid]]), 	
        z=t(ls_dat[[plotid]]),
        axes=FALSE,
        xlab="",
        ylab="",
        cex.lab=0.8,
        col=pal, 
        breaks=breaks, 
        main="")
    #mtext(text1, side=2, line=0.5, cex=0.7, font=2) # text vertical
    mtext(text1, side=2, line=0.5, cex=0.7, font=2,las=2) # text holizontal
    abline(v=c(12,53)+0.5,col="white",lwd=1,xpd=FALSE) # lines that devide sample groups
  }


  png(file="heatamap_MLL2_7.png",width = 8, height = 8, units="in",res=300)
  layout(matrix(data=c(0,1,2,3,0,4,5,0,6,7,0,8,9), nrow=13, byrow=FALSE), widths=c(1), heights=c(0.4,2,2,2,rep(c(0.4,2,2),3)), respect=FALSE)
  for (i in 1:length(textlab)){	
    meth_heatmap_2(ls_dat=dat.split.sorted, plotid=i, text1=textlab[i])
    cat(paste(i,"\n"))}
  dev.off()



  # Wilcoxon-ranksum test between wild tpye and mutant for each category
  # perhaps this really should be done via mt.teststat!
  Wilcox.test <- function(x,s1,s2) {
    x1 <- x[s1]
    x2 <- x[s2]
    x1 <- as.numeric(x1)
    x2 <- as.numeric(x2)
    wx.out <- wilcox.test(x1,x2, alternative="two.sided", paired = FALSE)
    p.Wilcox <- as.numeric(wx.out$p.value)
    return(p.Wilcox)
  }

  require(qvalue)
  qvalues_ls <- vector(mode="list", length=length(dat.split))
  names(qvalues_ls) <- names(dat.split)
  for (i in 1:length(dat.split)) {
    cat("\nPROBE CATEGORY:\n")
    print(names(dat.split)[[i]])
    rawp.wx <- apply(dat.split[[i]],1,Wilcox.test, s1=MutID, s2=WtID)
    q.mut <- qvalue(rawp.wx)
    qvalues_ls[[i]] <- q.mut
    summary(qvalues_ls[[i]])}
  # 
  # PROBE CATEGORY:
  # [1] "cgi_p"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	1	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value      7     50   448   1045  1855 3671 61309
  # q-value      0      0     0      0     0    0     7
  # 
  # 
  # PROBE CATEGORY:
  # [1] "cgi_p.PRC"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.8697504	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value     15     93   815   1839  3020 5144 31658
  # q-value      0      0     0      0     0    0 31880
  # 
  # 
  # PROBE CATEGORY:
  # [1] "non.cgi_p"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.742423	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05  <0.1    <1
  # p-value     86    714  4520   8960 14025 21629 93078
  # q-value      0      0     0      0     0   890 93615
  # 
  # 
  # PROBE CATEGORY:
  # [1] "cgi_ex"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.5712687	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value     44    205   943   1583  2206 3154 10340
  # q-value      0      2    17    162   665 1683 10383
  # 
  # 
  # PROBE CATEGORY:
  # [1] "non.cgi_ex"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.657696	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value     74    441  2265   3958  5756 8296 27361
  # q-value      0      0     0     74   885 3186 27515
  # 
  # 
  # PROBE CATEGORY:
  # [1] "cgi_ir"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.5851767	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value     35    221   998   1728  2497 3533 12137
  # q-value      0      0     7    111   553 1632 12187
  # 
  # 
  # PROBE CATEGORY:
  # [1] "non.cgi_ir"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.6410972	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05  <0.1    <1
  # p-value    208   1270  6988  12494 18311 26516 87699
  # q-value      0      0     0    165  2182 10023 88137
  # 
  # 
  # PROBE CATEGORY:
  # [1] "cgi_inter"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.6313061	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
  # p-value     65    281  1363   2432  3588 5148 17742
  # q-value      0      1     4    108   446 1879 17819
  # 
  # 
  # PROBE CATEGORY:
  # [1] "non.cgi_inter"
  # 
  # Call:
  # qvalue(p = rawp.wx)
  # 
  # pi0:	0.4734747	
  # 
  # Cumulative number of significant calls:
  # 
  #         <1e-04 <0.001 <0.01 <0.025 <0.05  <0.1    <1
  # p-value    248   1714  8774  15540 22069 30378 78034
  # q-value      0      0     1   2866 12181 26159 78332
  # 

}


