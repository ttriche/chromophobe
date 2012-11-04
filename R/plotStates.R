## plot emission and transition probabilities by state
##
plotStates <- function(object, what='emissions', states=NULL, target=NULL,...){
  if(what == 'emissions') {
    plotEmissions(emissions(object), states=states, target=target, ...)
  } else if(what == 'transitions') {
    plotTransitions(transitions(object), states=states, ...)
  } else { 
    stop(paste("Don't know how to plot", what))
  }
}

plotEmissions <- function(emissions, states=NULL, target=NULL, ...) { # {{{

  require('reshape2')
  require('pheatmap')

  ## deal with state assignments here:
  if(!is.null(states)) message('Have not dealt with state assignments yet :-/')

  ## now cast to a manageable structure and plot
  mat <- acast(emissions, state ~ mark, value.var='p')
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
  return(p)

} # }}}

plotTransitions <- function(trans, states=NULL, cluster=F, ...) { # {{{
  message("FIXME: Include mark frequencies as row/column annotations!")
  require(pheatmap)
  if(!is.null(states)) stop('Collapsing states is not yet supported')
  pheatmap(t(trans), cluster_col=cluster, cluster_row=cluster,
           ## add annotations here
           main='Transition probabilities, from state y to state x')
} # }}}
