setMethod("sort", signature(x="SummarizedExperiment"), 
          function(x) x[ names(sort(rowData(x))), ] ) 
