#' Perform tSNE
#'
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta simulation parameters
#' @param data expression matrix
#' @param plotname the name of the jpeg file
#' @param label the column name of the meta data that the points needs to be colored by
#' @import Rtsne
#' @import ggplot2
#' @export
PlotTsne <- function(meta, data, cif_type, pca = T, n_pc, perplexity=30, label, saving=F, plotname,system.color=T){
  uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[,uniqcols];meta <- meta[uniqcols,,drop=FALSE]
  uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
  data <- data[uniqrows,]
  if (pca){
    data_pc <- prcomp(t(data))
    data <- t(data_pc$x[,c(1:n_pc)])
  }
  data_tsne=Rtsne(t(data),perplexity=perplexity)
  plot_tsne <- cbind(meta, data.frame(label=factor(meta[,label]),x=data_tsne$Y[,1],y=data_tsne$Y[,2]))
  p <- ggplot(plot_tsne, aes(x, y))
  p <- p + geom_point(aes(colour = label),shape=20) + labs(color=label)
  p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

  if (system.color==F){
    if(cif_type=="discrete" | cif_type=="one.population"){
      color_5pop <- c("#F67670", "#0101FF", "#005826", "#A646F5", "#980000")
      names(color_5pop) <- 1:5
    }else{
      color_5pop <- c("#CC9521", "#1EBDC8", "#0101FF", "#005826", "#7CCC1F", "#A646F5", "#980000", "#F67670")
      names(color_5pop) <- c("6_7","7_8","8_2","8_3","7_9","9_4","9_5","6_1")
    }
    p <- p + scale_color_manual(values=color_5pop[levels(plot_tsne[['label']])])
  }
  if(saving==T){ggsave(p,filename=plotname,device='pdf',width=5,height=4)}
  if(saving==F){p <- p + ggtitle(plotname)}
  return(list(plot_tsne,p))
}



#' Perform UMAP
#'
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs UMAP, saves the result to a jpeg file with user specified name
#' @param meta simulation parameters
#' @param data expression matrix
#' @param plotname the name of the jpeg file
#' @param min_dist minimum distance to identify neighborhood,default is 0.5
#' @param n_neighbors The size of local neighborhood for manifold learning, default is 30
#' @param label the column name of the meta data that the points needs to be colored by
#' @import umap
#' @import ggplot2
#' @export
PlotUmap <- function(meta, data, pca = T, n_pc, min_dist=0.5, n_neighbors = 30,label, saving=F, plotname,system.color=T){
  custom.config = umap.defaults
  custom.config$min_dist = min_dist
  custom.config$n_neighbors = n_neighbors

  uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[,uniqcols];meta <- meta[uniqcols,,drop=FALSE]
  uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
  data <- data[uniqrows,]
  if (pca){
    data_pc <- prcomp(t(data),scale. = TRUE)
    data <- t(data_pc$x[,c(1:n_pc)])
  }
  umap_res <-umap(t(data),custom.config)
  umap_data <- umap_res$layout
  colnames(umap_data) <- c("x", "y")
  plot_umap <- cbind(meta, data.frame(label=factor(meta[,label])),umap_data)
  p <- ggplot(plot_umap, aes(x, y))
  p <- p + geom_point(aes(colour = label),shape=20) + labs(color=label)
  p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  if(saving==T){ggsave(p,filename=plotname,device='pdf',width=5,height=4)}
  if(saving==F){p <- p + ggtitle(plotname)}
  return(list(plot_umap,p))
}
