fromChr <- function(seqs, prefix='chr') { # {{{
  for(i in rev(seq_len(nchar(prefix)))) {
    seqs <- gsub(paste0('^', substr(prefix, 1, i)), '', seqs)
  }
  return(seqs)
} # }}}
toChr <- function(seqs, prefix='chr') { # {{{
  paste0(prefix, fromChr(seqs, prefix))
} # }}}
df2GR <- function(df, keepColumns=F, ignoreStrand=F, prefix='chr') { ## {{{

  require(GenomicRanges)
  if(class(df) == 'DataFrame') df <- as(df, 'data.frame')
  if(class(df) != 'data.frame') stop('df must be a data.frame or DataFrame')

  ## tidy up column names to coerce
  subs <- c(chr='chrom',
            seqnames='chrom',
            start='chromStart', 
            end='chromEnd')
  
  if(!all(unique(subs) %in% names(df))) {
    for(s in names(subs)) {
      names(df) = gsub(paste0('^', s, '$'), subs[s], names(df), ignore=TRUE)
    }
    if(!all(unique(subs) %in% names(df))) {
      stop('df must have columns chrom, chromStart, chromEnd to proceed')
    }
  }

  ## assign genome pre-emptively if possible
  if('genome' %in% names(attributes(df))) {
    g <- attr(df, 'genome') 
  } else {
    g <- NULL
  }

  ## fix seqnames if necessary 
  if(!all(grepl(prefix, df$chrom))) {
    df$chrom <- toChr(df$chrom, prefix=prefix)
  }

  ## fix starts/ends if necessary 
  df$chromStart <- as.numeric(df$chromStart)
  df$chromEnd <- as.numeric(df$chromEnd)

  ## fuss about any missing data, which will be dropped presently 
  if(any(is.na(df$chromStart))||any(is.na(df$chromEnd))||any(is.na(df$chrom))) {
    warning('Dropping ranges w/chrom/chromStart/chromEnd == NA')
    df <- df[ !is.na(df$chrom), ]
    df <- df[ !is.na(df$chromStart), ]
    df <- df[ !is.na(df$chromEnd), ]
  }

  ## default is NOT to ignore stranding, but we will if asked
  if(ignoreStrand == FALSE && ("strand" %in% names(df))) {
    if(is.numeric(df$strand)) {
      df$strand <- strandMe(df$strand)
    }
    GR <- with(df, GRanges(chrom, 
                           IRanges(start=as.numeric(chromStart), 
                                   end=as.numeric(chromEnd)), 
                           strand=strand))
  } else {
    GR <- with(df, GRanges(chrom, 
                           IRanges(start=as.numeric(chromStart), 
                                   end=as.numeric(chromEnd))))
  }

  ## were range names provided?
  if('name' %in% names(df)) {
    names(GR) <- df$name
    df$name <- NULL
  } else {
    names(GR) <- rownames(df)
  }

  ## keep metadata?
  if(keepColumns) {
    skipped = c("rangename","chrom","chromStart","chromEnd","width","strand")
    elementMetadata(GR) <- as(df[, setdiff(names(df), skipped), drop=F], 
                              "DataFrame")
  }

  ## chintzy hack from the original methLab
  if('X' %in% names(values(GR))) {
    if(all(is.na(GR$X))) {
      GR$X <- NULL
    } else {
      names(values(GR))[which(names(values(GR))=='X')]='score'
    }
  }

  ## assign genome to GR if known 
  if(!is.null(g)) genome(GR) <- g

  return(GR)
} # }}}
fetchGEOMethylationDataset <- function(GSE) { # {{{

  require(minfi)
  require(Biobase)
  require(GEOquery)
  gset <- getGEO(GSE)
  idx <- length(gset)
  platform <- annotation(gset[[1]])

  ## fixme: call combine() for multiple-gset same-platform GSEs
  if(length(gset) > 1) { # {{{ get the methylation component(s)
    idx <- grep("GPL13534", attr(gset, "names")) ## 450k
    if(!is.null(idx)) platform <- 'GPL13534'
    if(is.null(idx)) {
      idx <- grep("GPL8490", attr(gset, "names"))
      if(!is.null(idx)) platform <- 'GPL8490'
    }
  } # }}}

  gset <- gset[[idx]]
  sampleNames(gset) <- gset$title
  chrCol <- grep('^chr$', ignore.case=TRUE, names(fData(gset)))[1]
  startCol <- grep('^mapinfo$', ignore.case=TRUE, names(fData(gset)))[1]
  fData(gset)$chrom <- fData(gset)[ , chrCol]
  fData(gset)$chromStart <- fData(gset)[ , startCol]
  fData(gset)$chromEnd <- fData(gset)[ , startCol]
  fData(gset)$strand <- '*'

  ## hm27 is annotated against hg18
  if(platform == 'GPL8490') { # {{{
    row.dat <- df2GR(fData(gset), keep=TRUE)
    genome(row.dat) <- 'hg18' 
  } # }}}

  ## hm450 is against hg19
  ## add simple SNP probe fix
  if(platform == 'GPL13534') { # {{{
    rsProbes <- grep('^rs', featureNames(gset), value=T)
    if(length(rsProbes) > 0) {
      ## the following were merged at some point... drop them 
      rsProbes <- rsProbes[ !rsProbes %in% c('rs13369115','rs5936512') ]
      ## rsLocs <- SNPlocs.Hsapiens.dbSNP.20120608:::rsid2loc(rsProbes)
      ## snpLocs.450k <- data.frame(chrom=gsub('^ch','chr',names(rsLocs)), 
      ##                            chromStart=rsLocs, chromEnd=rsLocs, 
      ##                            strand='*')
      ## rownames(snpLocs.450k) <- rsProbes
      ## save(snpLocs.450k, file="~/chromophobe/data/snpLocs.450k.rda")
      data(snpLocs.450k)
      fData(gset)[rsProbes, c('chrom','chromStart','chromEnd')] <- 
        snpLocs.450k[rsProbes, c('chrom','chromStart','chromEnd')]
    }
    row.dat <- df2GR(fData(gset), keep=TRUE)
    genome(row.dat) <- 'hg19' 
    ## length(row.dat)
    ## [1] 485575
  } # }}}

  row.dat <- keepSeqlevels(row.dat, paste0('chr', c(1:22, 'X', 'Y')))
  preprocessing <- c(rg.norm=paste0('See GEO ',GSE,' for details'),
                     minfi=NA, manifest=paste0('GEO ', platform))
  gset <- GenomicRatioSet(gr=row.dat, Beta=exprs(gset)[names(row.dat),],
                          pData=pData(gset), annotation=platform,
                          preprocessMethod=preprocessing)
  exptData(gset) <- SimpleList(GSE=GSE)
  return(gset)

} # }}}
