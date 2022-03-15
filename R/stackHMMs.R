#' build a poor man's stacked HMM
#' 
#' @param HMMs    a GRangesList of HMMs
#' @param djHMM   optional disjoint HMM (else one will be assembled)
#' @param cols    optional color key (else a default simpleCols will be used)
#' @param simple  simplify? (default is TRUE, you're on your own if it's FALSE)
#' @param only3   keep only `name`, `thick`, `itemRgb` as mcols? (FALSE)
#' @param verbose be verbose? (FALSE)
#' @param BPPARAM BiocParallel object (default is SerialParam()) 
#' 
#' @return        a GRanges of stacked HMM states (one or more per segment) 
#' 
#' @import        GenomicRanges
#' @import        BiocParallel
#' 
#' @export
stackHMMs <- function(HMMs, djHMM=NULL, cols=NULL, simple=TRUE, only3=FALSE, verbose=FALSE, BPPARAM=SerialParam()) {

  # no/wrong colors provided?
  if (!simple & is.null(cols)) {
    stop("If not simplified, provide a non-null `cols` argument!")
  } else if (is.null(cols)) {
    cols <- .getSimpleCols()
  }

  # no disjoint HMM? 
  if (is.null(djHMM)) {
    djHMM <- .getDjHMM(HMMs, only3=only3, verbose=verbose, BPPARAM=BPPARAM)
  } else { 
    stopifnot("name" %in% names(mcols(djHMM)))
    if (!all(c("thick", "itemRgb") %in% names(mcols(djHMM)))) {
      djHMM$thick <- ranges(djHMM)
      djHMM$itemRgb <- .toHex("0,0,0")
    }
  }

  # both pure *and* mixed states
  splt <- split(djHMM, djHMM$name)
  mc <- names(mcols(djHMM))
  columns <- c("name","thick","itemRgb")
  if (only3) {
    mcc <- columns 
  } else { 
    mcc <- c(columns, setdiff(names(mcols(djHMM)), columns))
  }
 
  # easy case 
  if (verbose) message("Cleaning up segments in pure states... ")
  pure <- grep("/", names(splt), invert=TRUE, value=TRUE)
  ps <- GRangesList(bplapply(splt[pure], .aggPure, cols=cols, BPPARAM=BPPARAM))
  ps <- aggregateStates(unname(sort(unlist(ps))), only3=only3, BPPARAM=BPPARAM)
  mcols(ps) <- mcols(ps)[, mcc]
  if (verbose) message("Done.")

  # hard case
  mixed <- grep("/", names(splt), value=TRUE)
  if (verbose) message(length(mixed), " mixed states found. Expanding...")
  ms <- GRangesList(bplapply(splt[mixed], .aggMix, cols=cols, BPPARAM=BPPARAM))
  ms <- unname(sort(unlist(ms)))
  mcols(ms) <- mcols(ms)[, mcc]
  if (verbose) message("Done.")

  # merge both 
  sort(c(ps, ms))

}


# helper
.getSimpleCols <- function() {

  colrs <- c(Accessible="255,255,0",
             Bivalent="128,48,160",
             Enhancer="255,195,77",
             Het_Rpt_Qui="255,255,255",
             Promoter="255,0,0",
             Repressed="128,128,128",
             Transcribed="0,128,0")
  simpleCols <- data.frame(MNEMONIC=names(colrs), RGB=colrs)
  rownames(simpleCols) <- simpleCols$MNEMONIC
  return(simpleCols)

}


# helper, although maybe it should be its own exported function too
.getDjHMM <- function(HMMs, only3=FALSE, verbose=FALSE, clarify=FALSE, BPPARAM=SerialParam()) {

  if (verbose) message("Disjoining HMM states... ", appendLF=FALSE)
  dj <- disjoin(unlist(GRangesList(HMMs)))
  if (verbose) message("Done.")
  for (i in names(HMMs)) {
    if (verbose) message("Adding states from ", i)
    ol <- findOverlaps(dj, HMMs[[i]])
    mcols(dj)[, i] <- as.character(mcols(HMMs[[i]])[subjectHits(ol), "name"])
  }
  contigs <- unique(seqnames(dj))
  dj <- sort(unlist(GRangesList(bplapply(split(dj, seqnames(dj))[contigs], 
                                         .combStates, 
                                         verbose=verbose,
                                         clarify=clarify,
                                         BPPARAM=BPPARAM)), use.names=FALSE))
  aggregateStates(dj, only3=only3, verbose=verbose, BPPARAM=BPPARAM)

}


# helper
.combStates <- function(x, verbose=FALSE, collapse="/", clarify=FALSE) {

  contig <- unique(seqnames(x))
  if (verbose) message("Labeling states on ", contig, "... ", appendLF=FALSE)
  x$name <- 
    apply(mcols(x), 1, .combineStates, collapse=collapse, clarify=clarify)
  if (verbose) message("Done.")
  return(x)

}


# helper; also shorten & clarify states if requested
.combineStates <- function(x, collapse="/", clarify=FALSE) {

  states <- sort(unique(x))
  if (clarify & length(states) > 1) {
    states <- substr(states, 1, 3)
    subs <- c(Pro="Tss", Rep="Prc", Tra="Txn")
    for (s in names(subs)) states <- sub(s, subs[s], states)
  }
  paste(states, collapse=collapse)

}


# helper
.aggPure <- function(gr, cols, verbose=FALSE) {

  p <- unique(gr$name)
  idx <- match(p, cols$MNEMONIC)
  if (verbose) message("Aggregating pure ", p, " states... ", appendLF=FALSE)
  mcols(gr)[, "itemRgb"] <- .toHex(as.character(cols[idx, "RGB"]))
  mcols(gr)[, "thick"] <- ranges(gr)
  if (verbose) message("Done.")
  return(gr)

}


# helper 
.aggMix <- function(gr, cols, verbose=FALSE) {

  ms <- unique(gr$name)
  states <- strsplit(ms, "/")[[1]]
  idx <- match(states, cols$MNEMONIC)
  hex <- .toHex(as.character(cols[idx, "RGB"]))
  names(hex) <- states
  origlen <- length(gr)
  ns <- length(states)
  
  grr <- rep(gr, ns)
  columns <- c("name","thick","itemRgb")
  mcols(grr)[, "name"] <- rep(states, each=origlen)
  mcols(grr)[, "thick"] <- ranges(grr)
  mcols(grr)[, "itemRgb"] <- hex[grr$name]
  
  if (verbose) message("Expanding ", ms, " states... ", appendLF=FALSE)
  sgrr <- split(grr, grr$name)
  for (state in names(sgrr)) {
    runs <- setdiff(names(mcols(gr)), columns)
    wr <- apply(mcols(gr)[, runs], 1, .collapseStates, state=state)
    mcols(sgrr[[state]])$name <- paste0(mcols(sgrr[[state]])$name, ":", wr)
  }
  if (verbose) message("Done.")
  unlist(sgrr, use.names=FALSE)

}


# helper
.collapseStates <- function(x, state) {

  paste(names(which(x == state)), collapse="/")

}
