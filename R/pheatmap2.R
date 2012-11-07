jet.colors <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
                                     "#E0F3F8", "#91BFDB", "#4575B4")))(100) 

pheatmap2 <- function(mat, color=jet.colors, kmeans_k=NA, breaks=NA, 
                      border_color="grey60", cellwidth=NA, cellheight=NA, 
                      scale="none", cluster_rows=TRUE, cluster_cols=TRUE, 
                      clustering_distance_rows="euclidean", 
                      clustering_distance_cols="euclidean", 
                      clustering_method="complete", 
                      treeheight_row=ifelse(cluster_rows, 50, 0), 
                      treeheight_col=ifelse(cluster_cols, 50, 0), 
                      legend=TRUE, legend_breaks=NA, legend_labels=NA, 
                      annotation=NA, annotation_colors=NA, annotation_legend=T,
                      drop_levels=TRUE, show_rownames=T, show_colnames=T, 
                      main=NA, fontsize=10, fontsize_row=fontsize, 
                      fontsize_col=fontsize, display_numbers=F, 
                      number_format="%.2f", fontsize_number=0.8 * fontsize, 
                      filename=NA, width=NA, height=NA, ...) { # {{{

  mat <- as.matrix(mat)
  mat <- scale_mat(mat, scale)
  if(!is.na(kmeans_k)) {
    km <- kmeans(mat, kmeans_k, iter.max=100)
    mat <- km$centers
    t <- table(km$cluster)
    rownames(mat) <- sprintf("cl%s_size_%d", names(t), t)
  } else { 
    km <- NA
  }

  if(cluster_rows) {
    tree_row <- cluster_mat(mat, 
                            distance=clustering_distance_rows, 
                            method=clustering_method)
    mat <- mat[tree_row$order, ]
  } else {
    tree_row=NA
    treeheight_row=0
  }

  if(cluster_cols) {
    tree_col <- cluster_mat(t(mat), 
                            distance=clustering_distance_cols, 
                            method=clustering_method)
    mat <- mat[, tree_col$order]
  } else {
    tree_col=NA
    treeheight_col=0
  }

  if(display_numbers) {
    fmat=matrix(sprintf(number_format, mat), nrow=nrow(mat), ncol=ncol(mat))
    attr(fmat, "draw")=TRUE
  } else {
    fmat=matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
    attr(fmat, "draw")=FALSE
  }

  if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])) {
    if(length(legend_breaks) != length(legend_labels)) {
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }

  if(is.na(breaks[1])) {
    breaks=generate_breaks(as.vector(mat), length(color))
  }

  if(legend & is.na(legend_breaks[1])) {
      legend=grid.pretty(range(as.vector(breaks)))
      names(legend)=legend
  } else if(legend & !is.na(legend_breaks[1])) {
    legend <- legend_breaks[legend_breaks>=min(breaks) &
                            legend_breaks<=max(breaks)]
    if(!is.na(legend_labels[1])) {
      legend_labels <- legend_labels[legend_breaks>=min(breaks) & 
                                     legend_breaks <= max(breaks)]
      names(legend) <- legend_labels
    } else {
      names(legend)=legend
    }
  } else {
    legend=NA
  }
  mat <- scale_colours(mat, col=color, breaks=breaks)

  if(!is.na(annotation[[1]][1])) {
      annotation <- annotation[colnames(mat), , drop=F]
      annotation_colors <- generate_annotation_colours(annotation, 
                                                       annotation_colors, 
                                                       drop=drop_levels)
  }

  if(!show_rownames) rownames(mat)=NULL
  if(!show_colnames) colnames(mat)=NULL
  
  pheatmap_motor(mat, border_color=border_color, cellwidth=cellwidth, 
                 cellheight=cellheight, treeheight_col=treeheight_col, 
                 treeheight_row=treeheight_row, tree_col=tree_col, 
                 tree_row=tree_row, filename=filename, width=width, 
                 height=height, breaks=breaks, color=color, legend=legend, 
                 annotation=annotation, annotation_colors=annotation_colors, 
                 annotation_legend=annotation_legend, fontsize=fontsize, 
                 fontsize_row=fontsize_row, fontsize_col=fontsize_col, 
                 fmat=fmat, fontsize_number=fontsize_number, main=main, ...)
  invisible(list(tree_row=tree_row, tree_col=tree_col, kmeans=km))

} # }}}

pheatmap.motor <- function(matrix, border_color, cellwidth, cellheight, 
                           tree_col, tree_row, treeheight_col, treeheight_row,
                           filename, width, height, breaks, color, legend, 
                           annotation, annotation_colors, annotation_legend, 
                           main, fontsize, fontsize_row, fontsize_col, fmat, 
                           fontsize_number, ...) { # {{{

  grid.newpage()
  mindim=lo(coln=colnames(matrix), rown=rownames(matrix), 
            nrow=nrow(matrix), ncol=ncol(matrix), cellwidth=cellwidth, 
            cellheight=cellheight, treeheight_col=treeheight_col, 
            treeheight_row=treeheight_row, legend=legend, annotation=annotation,
            annotation_colors=annotation_colors,
            annotation_legend=annotation_legend,
            main=main, fontsize=fontsize, fontsize_row=fontsize_row, 
            fontsize_col=fontsize_col, ...)

  if (!is.na(filename)) {
    pushViewport(vplayout(1:5, 1:5))
    if (is.na(height)) {
      height=convertHeight(unit(0:1, "npc"), "inches", valueOnly=T)[2]
    }
    if (is.na(width)) {
      width=convertWidth(unit(0:1, "npc"), "inches", valueOnly=T)[2]
    }
    r=regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) stop("Improper filename")
    ending=substr(filename, r + 1, r + attr(r, "match.length"))
    f <- switch(ending, 
                pdf=function(x, ...) pdf(x, ...), 
                png=function(x, ...) png(x, units="in", res=300, ...), 
                jpeg=function(x, ...) jpeg(x, units="in", res=300, ...), 
                jpg=function(x, ...) jpeg(x, units="in", res=300, ...), 
                tiff=function(x, ...) tiff(x, units="in",res=300, ...), 
                bmp=function(x, ...) bmp(x, units="in", res=300, ...), 
                stop("File type should be: pdf, png, bmp, jpg, tiff"))
    f(filename, height=height, width=width)
    pheatmap_motor(matrix, cellwidth=cellwidth, cellheight=cellheight, 
                   border_color=border_color, tree_col=tree_col, 
                   tree_row=tree_row, treeheight_col=treeheight_col, 
                   treeheight_row=treeheight_row, breaks=breaks, 
                   color=color, legend=legend, annotation=annotation, 
                   annotation_colors=annotation_colors, 
                   annotation_legend=annotation_legend, 
                   filename=NA, main=main, fontsize=fontsize, 
                   fontsize_row=fontsize_row, fontsize_col=fontsize_col, 
                   fmat=fmat, fontsize_number=fontsize_number, ...)
    dev.off()
    upViewport()
    return()
  }

  if (mindim < 3) border_color=NA

  if (!is.na(main)) {
    pushViewport(vplayout(1, 2))
    draw_main(main, fontsize=1.3 * fontsize, ...)
    upViewport()
  }

  if (!is.na(tree_col[[1]][1]) & treeheight_col != 0) {
    pushViewport(vplayout(2, 2))
    draw_dendrogram(tree_col, horizontal=T)
    upViewport()
  }

  if (!is.na(tree_row[[1]][1]) & treeheight_row != 0) {
    pushViewport(vplayout(4, 1))
    draw_dendrogram(tree_row, horizontal=F)
    upViewport()
  }

  pushViewport(vplayout(4, 2))
  draw_matrix(matrix, border_color, fmat, fontsize_number)
  upViewport()

  if (length(colnames(matrix)) != 0) {
    pushViewport(vplayout(5, 2))
    pars=list(colnames(matrix), fontsize=fontsize_col, ...)
    do.call(draw_colnames, pars)
    upViewport()
  }

  if (length(rownames(matrix)) != 0) {
    pushViewport(vplayout(4, 3))
    pars=list(rownames(matrix), fontsize=fontsize_row, 
        ...)
    do.call(draw_rownames, pars)
    upViewport()
  }

  ## Hui's patch
  draw_annotations <- function(converted_annotations, border_color) {
    n=ncol(converted_annotations)
    m=nrow(converted_annotations)
    x=(1:m)/m - 1/2/m
    y=cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
    for (i in 1:m) {
      grid.rect(x=x[i], 
                unit(y[1:n], "bigpts"), 
                width=1/m, height=unit(8, "bigpts"), 
                gp=gpar(fill=converted_annotations[i,], col=border_color))
    }
    grid.text(names(annotation), x=1, unit(y[1:n], "bigpts"), just="right")
  }

  if (!is.na(annotation[[1]][1])) {
    pushViewport(vplayout(3, 2))
    converted_annotation=convert_annotations(annotation, annotation_colors)
    draw_annotations(converted_annotation, border_color)
    upViewport()
  }

  if (!is.na(annotation[[1]][1]) & annotation_legend) {
    if (length(rownames(matrix)) != 0) {
      pushViewport(vplayout(4:5, 5))
    } else {
      pushViewport(vplayout(3:5, 5))
    }
    draw_annotation_legend(annotation, annotation_colors, border_color, 
                           fontsize=fontsize, ...)
    upViewport()
  }

  if (!is.na(legend[1])) {
    length(colnames(matrix))
    if (length(rownames(matrix)) != 0) {
      pushViewport(vplayout(4:5, 4))
    } else {
      pushViewport(vplayout(3:5, 4))
    }
    draw_legend(color, breaks, legend, fontsize=fontsize, ...)
    upViewport()
  }
} # }}}
