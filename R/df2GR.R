require(GenomicRanges)
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
setAs("data.frame", "GRanges", function(from) df2GR(from, ...))
setAs("DataFrame", "GRanges", function(from) df2GR(from, ...))
as.data.frame.DataFrame <- selectMethod("as.data.frame", "DataFrame")
