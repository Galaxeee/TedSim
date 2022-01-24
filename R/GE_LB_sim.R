#' Sample Cell Identity Factors via Brownian motion, based on the cell lineage
#' @param par current parent node label
#' @param depth current depth of the parent node
#' @param anc_state Cell Identity Factor value of the parent node
#' @param edges edge list of the cell lineage tree
#' @param muts the mutation barcode matrix
#' @param p_d whether or not to have dropouts
#' @param mu mutation rate
#' @param edges_state edge list of the state tree
#' @param sif_mean state identity factor values
#' @param S state and state identity factor of the parent
#' @param p_a asymmetric division rate
#' @param cif cell identity factor matrix
#' @param flag whether or not generate barcode along with gene expressions
#' @param barcode current barcode of the parent node
#' @param N_ms Number of possible mutated states
#' @param unif_on if unif_on is TRUE, the mutated states will be synthetically generated using uniform distribution; otherwise, it will be sampled from a real dataset
#' @param lambda a num vector that indicates the value of lambda that weights the additional random walk value for different depths
Samplelineage <- function(par,depth,anc_state,edges, muts = NULL,p_d = 0.1, mu = 0.1,edges_state,sif_mean=NULL,S = NULL, p_a = 0.8,cif=NULL,flag=NULL, barcode = NULL, N_ms = NULL, unif_on = FALSE, lambda = NULL){
  children <- edges[edges[,2]==par,3] # get the children of the current node
  result<-lapply(c(1:length(children)),function(j){
    edge<-edges[edges[,2]==par & edges[,3]==children[j],] # given the parent and child, find the edge
    if(sum(edges[,2]==children[j])==0){
      result <- SampleEdgeNew(edge,depth,anc_state,edges,sif_mean = sif_mean,S=S,cif = cif, mu = mu,p_d = p_d, barcode = barcode, flag = flag,N_ms = N_ms, unif_on = unif_on,lambda = lambda)
      result <- result[c(1:(length(result[,1]-1))),]
    }else{
      result <- SampleEdgeNew(edge,depth,anc_state,edges,sif_mean = sif_mean,S=S,cif = cif, mu = mu,p_d = p_d, barcode = barcode, flag = flag,N_ms = N_ms, unif_on = unif_on,lambda = lambda)
      anc_state <- result[length(result[,1]),4]
      if (flag ==1){
        barcode <- result[length(result[,1]),5:length(result[1,])]
      }else{
        barcode <- NULL
      }
      result <- result[c(1:(length(result[,1]-1))),]

      depth <- depth + edge[4]
      result1 <- Samplelineage(children[j],depth,anc_state,edges,edges_state=edges_state,sif_mean=sif_mean,S = S, p_a = p_a,cif=cif,flag = flag,mu = mu, p_d = p_d,barcode = barcode,N_ms = N_ms, unif_on = unif_on,lambda = lambda)
      result <- rbind(result,result1)
    }
    return(result)
  })
  result<-do.call(rbind,result)
  result <-result[!duplicated(result[,4]),]
  result<- result[order(result[,2]),]
  return(result)
}

#' Generate child CIF based on parent CIF
#' @param edge current edge on the cell lineage tree
#' @param depth current depth of the parent node
#' @param anc_state Cell Identity Factor value of the parent node
#' @param edges edge list of the cell lineage tree
#' @param sif_mean state identity factor values
#' @param S state and state identity factor of the parent
#' @param cif cell identity factor matrix
#' @param mu mutation rate
#' @param p_d whether or not to have dropouts
#' @param barcode current barcode of the parent node
#' @param flag whether or not generate barcode along with gene expressions
#' @param N_ms Number of possible mutated states
#' @param unif_on if unif_on is TRUE, the mutated states will be synthetically generated using uniform distribution; otherwise, it will be sampled from a real dataset
#' @param lambda a num vector that indicates the value of lambda that weights the additional random walk value for different depths
SampleEdgeNew <- function(edge,depth,anc_state,edges,sif_mean=NULL,S=NULL,cif = NULL, mu = 0.1,p_d = 0, barcode = NULL, flag = 0, N_ms = NULL, unif_on = FALSE, lambda = NULL){
  state_prev <- S[edge[2],]
  state <- S[edge[3],]
  cifs <- sif_mean[[cif]]
  state_mean <-cifs[(cifs[,1]==state[1])&(cifs[,2]==state[2])&(cifs[,3]==state[3]),4]
  state_mean_prev <- cifs[(cifs[,1]==state_prev[1])&(cifs[,2]==state_prev[2])&(cifs[,3]==state_prev[3]),4]
  t_sample<-c(0,seq(0, edge[4], edge[4]))
  t_interval<-diff(t_sample)
  x_change <- sapply(t_interval,function(sig){rnorm(1,0,sqrt(sig))})[2]
  x_sample <- state_mean-state_mean_prev+lambda[depth+1]*cumsum(x_change)
  if (flag ==1){
    child_barcode <- barcode
    state_dist <- Mutated_state_dist(N_ms, cm)
    child_barcode <- generate_mutation(child_barcode,mu = mu,p_d = p_d,N_ms = N_ms, mutation_dist = state_dist, unif_on = unif_on)
    barcodes <- rbind(barcode,child_barcode)
    result<-cbind(depth+t_sample[-1],anc_state+x_sample,barcodes)
  }else{
    result<-cbind(depth+t_sample[-1],anc_state+x_sample)
  }
  rownames(result) <- c()
  result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
  return(result)
}

#' Generate joint profile of gene expression and lineage barcode for output
#' @param observed_counts Observed gene expression matrix
#' @param muts lineage character matrix
#' @param states cell states
#' @import DrImpute
#' @export
Generate_profile_multi <- function(observed_counts, muts, states){
  allele <- c()
  m_l <- c('D','I')
  for (k in 1:40){
    mutator <- m_l[sample(c(1,2),1)]
    num1 <- sample(1:100,1)
    num2 <- sample((2*num1):300,1)
    allele_temp <- sprintf("%g%s+%g",
                           num1, mutator,num2)
    allele <- c(allele,allele_temp)
  }

  counts <- observed_counts[[1]]
  exdata <- preprocessSC(counts)
  sf <- apply(exdata, 2, mean)
  npX <- t(t(exdata) / sf )
  lnpX <- log(npX+1)
  lnpX_imp <- DrImpute(lnpX)

  counts_t <- t(lnpX_imp)
  colnames(counts_t) <- paste("gene", seq(1, ncol(counts_t)), sep = "_")
  profiles <- data.frame(cell_id=character(), ClusterIdent = numeric(), HMID = character(),
                         stringsAsFactors = F)
  for (k in 1:ncells){
    barcode <- muts[k,]
    barcode <- replace(barcode, barcode=='-', '100D+100')
    barcode <- replace(barcode, barcode==0, 'None')
    for (j in 1:20){
      barcode <- replace(barcode, barcode==j, allele[j])
    }
    barcode <- paste(barcode,collapse = "-")
    cell_id <- paste("cell",states[k,4],sep = "_")
    #cell_id <- paste("f6_DEW",states_leaves[i,4],"bcAAVQ",sep = "_")
    temp_df <- data.frame(cell_id=cell_id,ClusterIdent = states[k,2], HMID = barcode,
                          stringsAsFactors = F)
    profiles <- rbind(profiles,temp_df)
  }

  profiles_out <- cbind(profiles,counts_t)

  var_gene <- c()
  for (gene in colnames(counts_t)){
    var_gene <- c(var_gene,var(counts_t[,gene]))
  }
  rank_gene <- data.frame(name = colnames(counts_t), var = var_gene)
  rank_gene <- rank_gene[order(-rank_gene$var),]
  colnames(rank_gene)<- c()
  returnList <- list(profiles_out,rank_gene)

  return(returnList)
}
