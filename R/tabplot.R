setMethod("tabplot", signature(dat="ANY"), # {{{
          function(dat, ...) {
            require(tabplot)
            tabplot::tableplot(dat=dat, ...)
          }) # }}}
setMethod("tabplot", signature(dat="data.frame"), # {{{
          function(dat, ...) {
            require(tabplot)
            tabplot::tableplot(dat=dat, ...)
          }) # }}}
setMethod("tabplot", signature(dat="data.table"), # {{{
          function(dat, ...) {
            require(tabplot)
            tabplot::tableplot(dat=dat, ...)
          }) # }}}
setMethod("tabplot", signature(dat="DataFrame"), # {{{
          function(dat, ...) {
            require(tabplot)
            tabplot::tableplot(dat=as.data.frame(dat), ...)
          }) # }}}
setMethod("tabplot", signature(dat="GRanges"), # {{{
          function(dat, ...) {
            require(tabplot)
            kept <- names(mcols(dat))
            if('dropcols' %in% names(list(...))) 
              kept <- with(list(...), setdiff(kept, dropcols))
            if('keepcols' %in% names(list(...))) 
              kept <- with(list(...), intersect(keepcols, kept))
            tabplot::tableplot(dat=as.data.frame(mcols(dat)[,kept]), ...)
          }) # }}}
