library(farver)

hex2rgb <- function(x) {
  apply(decode_colour(x), 1, paste, collapse=",")
}


rgb2hex <- function(x) {
  encode_colour(do.call(rbind, lapply(strsplit(decoded, ","), as.integer)))
}


stacked <- read.csv("stackHMM_state_annotations_processed.csv", row=1)
stacked$DESCRIPTION <- gsub("; +", "", stacked$DESCRIPTION)
stacked$COLOR <- toupper(stacked$COLOR)
stacked$RGB <- hex2rgb(stacked$COLOR)
simpleRGB <- c(Transcribed="0,128,0",
               Het_Rpt_Qui="255,255,255",
               Enhancer="255,195,77",
               Promoter="255,0,0",
               Bivalent="112,48,160",
               Repressed="128,128,128",
               CTCF="138,145,208")
stacked$RGBSIMPLE <- simpleRGB[stacked$SIMPLE]
stacked$STATE <- seq_len(nrow(stacked))
columns <- 
  c("STATE", "MNEMONIC", "RGB", "DESCRIPTION", "SIMPLE", "RGBSIMPLE", "NOTES")
stacked100state <- stacked[, columns]
write.table(stacked100state, "colors100state.tsv", 
            quote=FALSE, row.names=FALSE, sep="\t")

save(stacked100state, file="../../data/stacked100state.rda", compress="xz")
