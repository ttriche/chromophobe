if( 'regulatoR' %in% rownames(installed.packages()) ) {

  df2GR <- regulatoR:::df2GR

} else {

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
    for(s in names(subs)) names(df) = gsub(s, subs[s], names(df), ignore=TRUE)
    if(!all(unique(subs) %in% names(df))) {
      stop('df must have columns chrom, chromStart, chromEnd to proceed')
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

    ## fuss about any missing data, which will be dropped presently 
    if(any(is.na(df$chromStart))) warning('Dropping ranges w/chromStart == NA')
    if(any(is.na(df$chromEnd))) warning('Dropping ranges w/chromEnd == NA')
    if(any(is.na(df$chrom))) warning('Dropping ranges w/chrom == NA')
    df <- subset(df, !is.na(chromStart) & !is.na(chromEnd) & !is.na(chrom))

    ## default is NOT to ignore stranding, but we will if asked
    if(ignoreStrand == FALSE && ("strand" %in% names(df))) {
      if(is.numeric(df$strand)) {
        df$strand <- strandMe(df$strand)
      }
      GR <- with(df, GRanges(chrom, 
                             IRanges(start=chromStart, 
                                     end=chromEnd), 
                             strand=strand))
    } else {
      GR <- with(df, GRanges(chrom, 
                             IRanges(start=chromStart, 
                                     end=chromEnd)))
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

}

## add coercion
setAs("data.frame", "GRanges", function(from) df2GR(from, ...))
