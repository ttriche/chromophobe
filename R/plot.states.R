## plot emission and transition probabilities by state
##
plot.states <- function(model, what='emissions', states=NULL, target=NULL,...){
  require('reshape2')
  require('pheatmap')
  if(what == 'emissions') {
    if(!is.null(states)) {
      stop('Collapsing states need to be implemented')
    } else { 
      mat <- acast(model$emissions, state ~ mark, value.var='p')
    }
    if(!is.null(target)) {
      p <- pheatmap(mat, color=colorRampPalette(c('white','darkblue'))(100), 
                    legend=F, fontsize=16, cellwidth=16, cellheight=16,
                    clustering_distance_rows='manhattan', 
                    clustering_distance_cols='manhattan',
                    border_color='white', 
                    kmeans_k=target, ...)
    } else {
      p <- pheatmap(mat, color=colorRampPalette(c('white','darkblue'))(100), 
                    legend=F, fontsize=16, cellwidth=16, cellheight=16,
                    clustering_distance_rows='manhattan', 
                    clustering_distance_cols='manhattan',
                    border_color='white', ...)
    }
  } else {
    stop('Transition matrix plots are not yet supported')
  }
}
