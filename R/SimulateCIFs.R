#' Simulate Cell Identity Factor Matrix With Lineage Barcodes
#' @param ncells Number of Cells
#' @param N_states Number of leaf states to generate the cell state tree
#' @param cif_center Mean of CIFs, default is 1
#' @param Sigma Standard deviation of non-diff CIFs, difault is 0.5
#' @param p_a The asymmetric division rate, default is 0.8
#' @param p_edge The edge transition probability table, default is NULL so that the child edges are chosen with equal possibilities
#' @param n_CIF Total number of SIFs, both non-diff and diff combined
#' @param n_diff Number of diff SIFs
#' @param step Sampling stepsize of diff-SIF Brownian motion
#' @param p_d Dropout rate of CRISPR/Cas9 lineage barcodes
#' @param mu mutation rate of CRISPR/Cas9 lineage barcodes
#' @param N_char number of character sites for CRISPR/Cas9 lineage barcodes
#' @param N_ms number of mutated states for CRISPR/Cas9 lineage barcodes
#' @param unif_on sampling from synthetic uniform-distributed mutated states. When set FALSE, mutated states will be drawn from a real experimental dataset
#' @param SIF_res Optional input for State Identity Factors. If not, the function will generate SIF internally.
#' @param max_walk maximum walk distance of one asymmetric division on the cell state tree
#' @import ape
#' @export
SimulateCIFs <- function(ncells,phyla,cif_center=1, Sigma=0.5, p_a = 0.8,p_edge = NULL,n_CIF,n_diff,step = 1,p_d = 0.1,mu = 0.1, N_char = 9, N_ms = 100, unif_on = FALSE, SIF_res = NULL, max_walk = 2){
  T_cell <- stree(ncells,type = "balanced")
  T_cell$edge.length <- rep(1, length(T_cell$edge[,1]))
  N_nodes <- length(T_cell$edge[,1])+1
  cell_edges <- cbind(T_cell$edge,T_cell$edge.length)
  cell_edges <- cbind(c(1:length(cell_edges[,1])),cell_edges)
  cell_connections <- table(c(cell_edges[,2],cell_edges[,3]))
  cell_root <- as.numeric(names(cell_connections)[cell_connections==2])
  Node_cell <- cell_root

  if (is.null(SIF_res)){
    returnlist <- SIFGenerate(phyla,n_diff,step = step)
  }else{
    returnlist <- SIF_res
  }

  sif_mean <- returnlist$sif_mean
  sif_mean_raw <- returnlist$sif_mean_raw
  state_tree <- returnlist$tree
  state_edges <- cbind(state_tree$edge,state_tree$edge.length)
  state_edges <- cbind(c(1:length(state_edges[,1])),state_edges)
  state_connections <- table(c(state_edges[,2],state_edges[,3]))
  state_root <- as.numeric(names(state_connections)[state_connections==2])
  Node_state <- state_root

  sif_label <- sif_mean_raw[[1]][,2]

  S <- sif_mean[[1]]
  S <- S[S[,3]==0,]
  if (length(S)>4){
    S <- t(matrix(c(S[1,1:3],cell_root),byrow = TRUE))
  } else{
    S <- t(matrix(c(S[1:3],cell_root),byrow = TRUE))
  }

  State_table <- SimulateCellStates(cell_root,cell_edges,state_edges,sif_mean = sif_mean[[1]],S = S,p_a = p_a,p_edge,max_walk = max_walk)

  root_barcode <- rep(0,N_char)
  neutral <- Samplelineage(Node_cell,0,cif_center, edges = cell_edges, edges_state=state_edges,sif_mean=sif_mean,S = State_table,cif = 1, p_a = p_a,p_d = p_d, mu = mu, flag = 1, barcode = root_barcode, N_ms = N_ms, unif_on = unif_on)
  muts <- neutral[,5:length(neutral[1,])]

  param_names <- c("kon", "koff", "s")
  N_DE_cifs = c(0,0,n_diff)
  N_ND_cifs =c(n_CIF,n_CIF,n_CIF-n_diff)
  cifs <- lapply(c(1:3),function(parami){
    nd_cif <- lapply(c(1:N_ND_cifs[parami]),function(icif){
      rnorm(N_nodes-1,cif_center,Sigma)
    })
    nd_cif <- do.call(cbind,nd_cif)
    if(N_DE_cifs[parami]!=0){
      #if there is more than 1 de_cifs for the parameter we are looking at
      de_cif <- lapply(c(1:N_DE_cifs[parami]),function(cif_i){
        Samplelineage(Node_cell,0,cif_center,edges=cell_edges,edges_state=state_edges,sif_mean=sif_mean,S = State_table, p_a = p_a,cif=cif_i, flag = 0)
      })

      de_cif <- lapply(de_cif,function(X){X[,4]})
      de_cif <- do.call(cbind,de_cif)
      cifs <- cbind(nd_cif,de_cif)
      colnames(cifs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_cif[1,])),c(1:length(nd_cif[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_cif[1,])),c(1:length(de_cif[1,])),sep='_'))
    }else{
      cifs <- nd_cif
      colnames(cifs)<-paste(param_names[parami],rep('nonDE',length(nd_cif[1,])),c(1:length(nd_cif[1,])),sep='_')
    }
    return(cifs)
  })
  muts[muts == Inf] <- '-'
  colnames(State_table) <- c("parent","cluster","depth","cellID")
  cif_res <- list(cifs,State_table,state_tree,T_cell,sif_mean_raw,sif_label,muts)
  return(cif_res)

}
