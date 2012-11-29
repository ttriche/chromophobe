setMethod("keepSeqlevels", 
          signature(x="SummarizedExperiment", value="ANY"), 
          function(x, value) {
            y <- which(rownames(x) %in% names(keepSeqlevels(rowData(x),value)))
            return(x[y, ])
          })
