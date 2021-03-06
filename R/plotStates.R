## plot emission and transition probabilities by state
##
plotStates <- function(object, what='emissions', states=NULL, target=NULL,...){
  if(is.null(states)) states <- states(object)
  if(what == 'emissions') {
    plotEmissions(emissions(object), states=states, target=target, ...)
  } else if(what == 'transitions') {
    plotTransitions(transitions(object), states=states, ...)
  } else { 
    stop(paste("Don't know how to plot", what))
  }
}

reshapeEmissions <- function(emissions, how='matrix') { # {{{ refactored
  require('reshape2')
  if(how == 'matrix') {
    acast(emissions, state ~ mark, value.var='p')
  } else if(how == 'data.frame') {
    dcast(emissions, state ~ mark, value.var='p')
  } else {
    stop(paste("Don't know how to reshape emissions into a", how))
  }
} # }}}

plotEmissions <- function(emissions,states=NULL,target=NULL,cluster=T,legend=F, ...) { # {{{

  require('pheatmap')

  ## now cast to a manageable structure and plot
  if(!is.null(target)) {
    pheatmap2(reshapeEmissions(emissions, 'matrix'),
              color=colorRampPalette(c('white','darkblue'))(100), 
              legend=legend, fontsize=16, cellwidth=16, cellheight=16,
              clustering_distance_rows='manhattan', 
              clustering_distance_cols='manhattan',
              border_color='white', 
              kmeans_k=target, 
              main=paste('Emissions for', target, 'target states'),
              ... )
  } else {
    nStates <- length(unique(emissions$state))
    emis <- reshapeEmissions(emissions, 'matrix')[ states$Id, ] 
    rownames(emis) <- paste(states$Name, paste0('(', rownames(emis), ')'))
    pheatmap2(emis,
              color=colorRampPalette(c('white','darkblue'))(100), 
              legend=legend, fontsize=16, cellwidth=16, cellheight=16,
              cluster_rows=cluster, cluster_cols=cluster, 
              clustering_distance_rows='manhattan', 
              clustering_distance_cols='manhattan',
              main=paste('Emissions by state',ifelse(cluster,'(clustered)','')),
              border_color='white', ... )
  }

} # }}}

plotTransitions <- function(transitions, emissions=NULL, cluster=F, ...) { # {{{
  require(pheatmap)
  if(is.null(emissions)) {
    pheatmap2(t(transitions), cluster_col=cluster, cluster_row=cluster,
              fontsize=12, cellwidth=16, cellheight=16,
              main='Transition probabilities, from state y to state x', ...)
  } else { 
    ann_colors <- list()
    ann <- reshapeEmissions(emissions, 'data.frame')[,-1]
    for(i in names(ann)) ann_colors[[i]] <- c('white','darkblue')
    pheatmap2(t(transitions), cluster_col=cluster, cluster_row=cluster,
              fontsize=12, cellwidth=16, cellheight=16,
              annotation=ann, annotation_colors=ann_colors, annotation_legend=F,
              main=paste('Transition & emission probabilities by state'),
              ...)
  }
} # }}}
