#' Generate Cell Identity Factors recursively via Brownian motion
#' @param par current parent node label
#' @param depth current depth of the parent node
#' @param anc_state accumulated state identity factor value
#' @param edges edge list of the cell lineage
SampleCIF <- function(par,depth,anc_state,edges){
  resolution <- 0.01
  children <- edges[edges[,2]==par,3]
  result<-lapply(c(1:length(children)),function(j){
    edge<-edges[edges[,2]==par & edges[,3]==children[j],]
    if(sum(edges[,2]==children[j])==0){
      t_sample<-c(0,seq(0, edge[4], resolution))
      t_interval<-diff(t_sample)
      x_change <- sapply(t_interval,function(sig){rnorm(1,0,2 * sqrt(sig))})
      x_sample <- cumsum(x_change)
      x_sample_smooth <- smooth.spline(x_sample,df = 1/resolution)
      x_sample <- x_sample_smooth$y
      result<-cbind(depth+t_sample[-1],x_sample+anc_state)
      result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
      result <- result[c(1:(length(result[,1]-1))),]
    }else{
      t_sample<-c(0,seq(0, edge[4], resolution))
      t_interval<-diff(t_sample)
      x_change <- sapply(t_interval,function(sig){rnorm(1,0,2 * sqrt(sig))})
      x_sample <- cumsum(x_change)
      x_sample_smooth <- smooth.spline(x_sample,df = 1/resolution)
      x_sample <- x_sample_smooth$y
      result<-cbind(depth+t_sample[-1],x_sample+anc_state)
      result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
      result <- result[c(1:(length(result[,1]-1))),]
      anc_state <- result[length(result[,1]),4]
      result <- result[c(1:(length(result[,1]-1))),]
      depth <- depth + edge[4]
      result1 <- SampleCIF(children[j],depth,anc_state,edges)
      result <- rbind(result,result1)
    }
    return(result)
  })
  result<-do.call(rbind,result)
  result <-result[!duplicated(result[,4]),]
  result<- result[order(result[,2]),]
  return(result)
}


#' Simulate State Identity Factor values given a state tree
#' @param phyla state tree structure
#' @param n_diff number of diff-CIF simulated
#' @param step stepsize of sampling intermediate states
#' @import phytools
#' @import plyr
#' @import ape
#' @export
SIFGenerate <- function(phyla, n_diff, step = 1){
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
  connections <- table(c(edges[,2],edges[,3]))
  Node_par <- as.numeric(setdiff(levels(as.factor(edges[,2])), levels(as.factor(edges[,3]))))
  N_nodes <- length(phyla$edge[,1])+1
  M <- matrix(0,N_nodes,n_diff)
  de_sif <- lapply(c(1:n_diff),function(sif_i){
    SampleCIF(Node_par,0,1,edges)
  })
  de_sif_raw <- de_sif
  de_sif <- lapply(de_sif,function(X){
    X[,3] <- floor(round(X[,3]/step,digits=2)* 100) / 100
    X_sample <- as.data.frame(X[X[,3]%%1==0,])
    X_sample <- X_sample[!duplicated(X_sample[1:3]),]
    x_sample <- as.matrix(X_sample, rownames=TRUE)
    return(x_sample)
  })
  returnList <- list("tree" = phyla, "sif_mean" = de_sif, "sif_mean_raw" = de_sif_raw)
  return(returnList)
}


#' Find the child state based on parent state and step size
#' @param S_parent current state of the parent node
#' @param child child node label
#' @param state_edges edge list of the cell state tree
#' @param sif_mean edge list of the cell lineage
WalkTree <- function(S_parent,child,state_edges,sif_mean, p_edge, step){
  S_child <- c(S_parent[1:3],child)
  #state_par <- S_parent[1] #previous state of par
  #target_state <- S_parent[2] #target state of par
  #depth_par <- S_parent[3] #depth of par
  depth_range <- sif_mean[sif_mean[,1]==S_child[1]&sif_mean[,2]==S_child[2],3]
  depth_room <- depth_range[length(depth_range)] - S_child[3]

  while (step > 0){
    if (step >= depth_room){
      step <- step - depth_room
      child_states <- state_edges[state_edges[,2]==S_child[2],3]
      if(sum(state_edges[,2]==S_child[2])!=0){
         if (is.null(p_edge)){
           p_childs <- NULL
         }else{
           p_childs <- p_edge[p_edge[,1]==S_child[2],3]
         }

         S_child[1] <- S_child[2]
         choice <- sample(length(child_states),1,prob = p_childs)
         S_child[2] <- child_states[choice]
         depth_range <- sif_mean[sif_mean[,1]==S_child[1]&sif_mean[,2]==S_child[2],3]
         S_child[3] <- min(depth_range)
         depth_room <- depth_range[length(depth_range)] - S_child[3]
      }else{
        S_child[3] <- S_child[3] + depth_room
        step <- 0
      }
    }else{
      S_child[3] <- S_child[3] + step
      step <- 0
    }
  }
  return (S_child)
}

#' Simulate cell states based on the state tree
#' @param par current parent node
#' @param cell_edges edge list of the cell lineage tree
#' @param state_edges edge list of the cell state tree
#' @param sif_mean State Identity Factors
#' @param S current state of the parent node
#' @param p_a asymmetric division rate
#' @param p_edge branching possibilities
#' @param max_walk maximum walk distance on the state tree for one asymmetric division
SimulateCellStates <- function(par,cell_edges,state_edges,sif_mean,S,p_a,p_edge = NULL,max_walk = 2){
  children <- cell_edges[cell_edges[,2]==par,3] # get the children of the current node
  State_list <-S
  flag <-sample(c(1,2),1)
  State_list<-lapply(c(1:length(children)),function(j){
    #browser()
    edge<-cell_edges[cell_edges[,2]==par & cell_edges[,3]==children[j],] # given the parent and child, find the edge
    S_parent <- S

    state_par <- S_parent[1] #previous state of par
    target_state <- S_parent[2] #target state of par
    depth_par <- S_parent[3] #depth of par
    S_child <- c(S_parent[1:3],children[j])

    if (depth_par == 0){
      child_states <- state_edges[state_edges[,2]==state_par,3]
      if (is.null(p_edge)){
        p_childs <- NULL
      }else{
        p_childs <- p_edge[p_edge[,1]==state_par,3]
        if (sum(p_childs) != 1){
          print("error!")
        }
      }
      if (runif(1,0,1)<=p_a){
        if (j==flag){
          choice <- sample(length(child_states),1,prob = p_childs)
          S_parent[2] <-child_states[choice]
        }
      }
    }

    if (runif(1,0,1)<=p_a){
      if (j==flag){
        step_state <- sample(1:max_walk,1)
        S_child <- WalkTree(S_parent,children[j],state_edges,sif_mean, p_edge, step = step_state)
      }
    }


    #recursive calling if not leaf node
    if(sum(cell_edges[,2]==children[j])==0){
      #current node is leaf node
      State_list <- rbind(t(matrix(State_list)),t(matrix(S_child)))
    }
    else{
      result1 <- SimulateCellStates(children[j],cell_edges,state_edges,sif_mean = sif_mean,S = S_child,p_a = p_a,p_edge = p_edge, max_walk = max_walk)
      State_list <- rbind(t(matrix(State_list)),result1)
    }
    return(State_list)
  })
  State_list<-do.call(rbind,State_list)
  State_list <-State_list[!duplicated(State_list[,4]),]
  State_list<- State_list[order(State_list[,4]),]
  return(State_list)
}

#' Simulate true count matrix given simulated CIF values
#' @param ngenes Number of genes
#' @param ncif number of cifs simulated
#' @param ge_prob ge probability
#' @param ncells Number of cells simulated
#' @param cif_res the CIFs simulated
#' @param prob_hge the probability of hge, default is 0.015
#' @param mean_hge the mean of hge, default is 5
#' @param scale_s transcription rate scaler, or a vector to specify cell-type specific scale_s
#' @import phytools
#' @import plyr
#' @export
CIF2Truecounts <- function(ngenes,ncif,ge_prob,ncells, cif_res, prop_hge = 0.015, mean_hge = 5, scale_s = 1){
  seed <- sample(c(1:1e+05), size = 2)
  gene_effects <- GeneEffects(ngenes = ngenes, ncif = ncif,
                              randseed = seed[2], prob = ge_prob, geffect_mean = 0,
                              geffect_sd = 1)
  data(match_params)
  match_params[, 1] = log(base = 10, match_params[, 1])
  match_params[, 2] = log(base = 10, match_params[, 2])
  match_params[, 3] = log(base = 10, match_params[, 3])
  match_params_den <- lapply(c(1:3), function(i) {
    density(match_params[, i], n = 2000)
  })
  params <- Get_params(gene_effects, cif_res[[1]], match_params_den,
                       0, scale_s = scale_s)

  counts <- lapply(c(1:ngenes), function(i) {
    count <- sapply(c(1:ncells), function(j) {
      y <- rbeta(1, params[[1]][i, j], params[[2]][i,
                                                   j])
      x <- rpois(1, y * params[[3]][i, j])
      return(x)
    })
  })

  chosen_hge <- sample(ngenes, ceiling(ngenes * prop_hge),
                       replace = F)
  multi_factors <- numeric(length(chosen_hge))
  rank_sum <- rank(rowSums(params[[3]][chosen_hge, ]))
  multi_factors <- sapply(1:length(chosen_hge), function(igene) {
    tosubstract <- -rank_sum[igene] * 1/length(chosen_hge) + 1
    if (runif(1, 0, 1) < 1) {
      multi_factor <- mean_hge - tosubstract
    }
    else {
      multi_factor <- mean_hge
    }
    return(multi_factor)
  })
  new_s <- matrix(0, length(chosen_hge), ncells)
  for (i in 1:length(chosen_hge)) {
    new_s[i, ] <- params[[3]][chosen_hge[i], ] * (2^multi_factors[i])
  }
  params[[3]][chosen_hge, ] <- new_s
  counts[chosen_hge] <- lapply(1:length(chosen_hge), function(i) {
    s_vec <- new_s[i, ]
    count <- sapply(c(1:ncells), function(j) {
      x <- rpois(1, s_vec[j])
      return(x)
    })
    return(count)
  })
  chosen_hge <- cbind(chosen_hge, multi_factors)
  cell_meta <- cbind(cellid = paste("cell", seq(1, ncells),
                                    sep = "_"), cif_res[[2]][1:ncells,], cif_res[[1]])
  counts <- do.call(rbind, counts)

  true_counts_res <-list(counts = counts, gene_effects = gene_effects,
                         cell_meta = cell_meta, kinetic_params = params)
  return(true_counts_res)
}

#' Simulate observed count matrix given technical biases and the true counts
#' @param true_counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles, default is 16
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @param hge2true if we add high gene expression to true counts
#' @param SE input, should be a summerized experiment rather than a list of elements, default is False
#' @param nbatch number of batches
#' @import SummarizedExperiment
#' @export
True2ObservedCounts <- function(SE=NULL,true_counts,meta_cell,protocol,alpha_mean=0.1,alpha_sd=0.002,
                                lenslope=0.02,nbins=20,gene_len,amp_bias_limit=c(-0.2, 0.2),
                                rate_2PCR=0.8,nPCR1=16, nPCR2=10, LinearAmp=F, LinearAmp_coef=2000,
                                depth_mean, depth_sd, nbatch=1){
  if(!is.null(SE)){
    meta_cell <- colData(SE)
    true_counts <- assays(SE)$count
  }
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005; depth_lb <- 200 # lower bound for capture efficiency and sequencing depth
  rate_2cap_vec <- rnorm_trunc(n=ncells, mean = alpha_mean, sd=alpha_sd, a=rate_2cap_lb, b=Inf)
  depth_vec <- rnorm_trunc(n=ncells, mean = depth_mean, sd=depth_sd,a=depth_lb,b=Inf)

  observed_counts <- lapply(c(1:ncells),function(icell){
    amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol,
                  rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias,
                  rate_2PCR=rate_2PCR, nPCR1=nPCR1, nPCR2=nPCR2, LinearAmp = LinearAmp,
                  LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[icell])
  })
  ## assign random batch ID to cells
  batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
  meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec, batch=batchIDs)
  meta_cell <- cbind(meta_cell, meta_cell2)

  if (protocol=="UMI"){
    UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
    nreads_perUMI <- lapply(observed_counts, "[[", 2)
    nUMI2seq <- sapply(observed_counts, "[[", 3)
    observed_counts <- UMI_counts
  } else
    observed_counts <- do.call(cbind,observed_counts)

  ## add batch effects to observed counts
  # use different mean and same sd to generate the multiplicative factor for different gene in different batch
  if (nbatch>1){
    mean_matrix <- matrix(0, ngenes, nbatch)
    batch_effect_size <- 2
    gene_mean <- rnorm(ngenes, 0, 1)
    temp <- lapply(1:ngenes, function(igene) {
      return(runif(nbatch, min = gene_mean[igene]-batch_effect_size, max = gene_mean[igene]+batch_effect_size))
    })
    mean_matrix <- do.call(rbind, temp)

    batch_factor <- matrix(0, ngenes, ncells)
    for (igene in 1:ngenes){
      for (icell in 1:ncells){
        batch_factor[igene, icell] <- rnorm(n=1, mean=mean_matrix[igene, batchIDs[icell]], sd=0.01)
      }
    }
    observed_counts <- 2^(log2(observed_counts)+batch_factor)
  }

  if(is.null(SE)){
    if (protocol=="UMI"){return(list(counts=observed_counts, cell_meta=meta_cell, nreads_perUMI=nreads_perUMI,
                                     nUMI2seq=nUMI2seq))} else
                                       return(list(counts=observed_counts, cell_meta=meta_cell))
  } else{
    assays(SE)$observed_counts <- observed_counts
    colData(SE)<-meta_cell
    return(SE)
  }
}

#' Simulate technical biases
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit[2]/(nbins-1)
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
cal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit){

  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins))*lenslope
  len_bias_bin <- len_bias_bin-median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)

  rand_bias <- rnorm(ngenes, mean=0, sd=max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)

  binsize <- floor(ngenes/nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins-1)){
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin-1)*binsize+1) : (ibin*binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins-1)*binsize+1) : ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins

  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias+len_bias
  return(amp_bias)
}

#' expand transcript counts to a vector of binaries of the same length of as the number of transcripts
#' @param true_counts_1cell number of transcript in one cell
expand2binary <- function(true_counts_1cell){
  expanded_vec <- rep(1, sum(true_counts_1cell))
  trans_idx <- sapply(which(true_counts_1cell>0),
                      function(igene){return(rep(igene, true_counts_1cell[igene]))})
  trans_idx <- unlist(trans_idx)
  return(list(expanded_vec, trans_idx))
}

#' Getting GeneEffects matrices
#'
#' This function randomly generates the effect size of each cif on the dynamic expression parameters
#' @param ngenes number of genes
#' @param ncif number of cifs
#' @param randomseed (should produce same result if ngenes, ncif and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from
#' @return a list of 3 matrices, each of dimension ngenes * ncif
GeneEffects <- function(ngenes,ncif,randseed,prob,geffect_mean,geffect_sd){
  #set.seed(randseed)
  lapply(c('kon','koff','s'),function(param){
    effect <- lapply(c(1:ngenes),function(i){
      nonzero <- sample(size=ncif,x=c(0,1),prob=c((1-prob),prob),replace=T)
      nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind,effect))
  })
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from cif and gene effects
#'
#' This function takes gene_effect and cif, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has ncif columns, and ngenes rows.
#' @param cif a vector of length ncif, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples
#' Get_params()
Get_params <- function(gene_effects,cif,match_param_den,bimod,scale_s = 1){
  params <- lapply(1:3, function(iparam){cif[[iparam]] %*% t(gene_effects[[iparam]])})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2)
    # this line is to check that the row and columns did not flip
    temp <- alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })

  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x*scale_s}))

  return(scaled_params)
}


#' Getting the parameters for simulating gene expression from cif and gene effects
#'
#' This function takes gene_effect and cif, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param cif a vector of length ncif, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has ncif columns, and ngenes rows.
#' @param param_realdata the fitted parameter distribution to sample from
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples
#' Get_params()
Get_params2 <- function(gene_effects,cif,bimod,ranges){
  params <- lapply(gene_effects,function(X){cif %*% t(X)})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
}

#' @export
get_prob <- function(glength){
  if (glength >= 1000){prob <- 0.7} else{
    if (glength >= 100 & glength < 1000){prob <- 0.78}
    else if (glength < 100) {prob <- 0}
  }
  return(prob)
}

#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
#' @export
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias,
                          rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ){
  ngenes <- length(gene_len)
  if (protocol=="nonUMI"){data(len2nfrag)} else
    if(protocol=="UMI"){ } else
    {stop("protocol input should be nonUMI or UMI")}
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expand2binary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]

  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]

  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if (sum(captured_vec[1:(length(captured_vec)-1)]) < 1) {return(rep(0, ngenes))}
  captured_vec[length(captured_vec)] <- 1
  inds[[2]] <- which(captured_vec > 0); captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]

  amp_rate <- c((rate_2PCR+amp_bias[trans_idx[1:(length(trans_idx)-1)]]),1)

  # pre-amplification:
  if (LinearAmp){
    PCRed_vec <- captured_vec*LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp*2+captured_vec-temp
    for (iPCR in 2:nPCR1){
      eff <- runif(length(temp))*amp_rate
      v1 <- temp*(1-eff)
      round_down <- (v1-floor(v1)) < runif(length(v1))
      v1[round_down] <- floor(v1[round_down]); v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2*(temp-v1)
    }
    PCRed_vec <- temp
  }

  if (protocol=="nonUMI"){ # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    amp_mol_count=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]]
      amp_mol_count[i]=sum(x)
    }

    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count>0)){
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),],
                                    amp_mol_count[igene], replace = TRUE))}
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    for (iPCR in 3:nPCR2){
      frag_vec <- frag_vec + round(frag_vec*rate_2PCR)
    }
    SEQ_efficiency=N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1) {read_count <- frag_vec} else{
      read_count <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)}) }
    return(read_count)
  } else if (protocol=="UMI"){

    prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx)-1)]], get_prob)
    # fragmentation:
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], prob = prob_vec[igene] ))})

    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }

    frag_vec <- round(frag_vec * (1+rate_2PCR)^(nPCR2-1))

    SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
      sequenced_vec <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)})}

    temp_vec <- c(sequenced_vec,1)
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]

    UMI_counts=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]];
      UMI_counts[i]=sum(x>0);
    }

    return(list(UMI_counts, sequenced_vec, sum(frag_vec>0)))
  }
}

#' sample from truncated normal distribution
#' @param a the minimum value allowed
#' @param b the maximum value allowed
rnorm_trunc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a | temp > b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}


#' Simulate observed count matrix given technical biases and the true counts
#' @param phyla cell state tree
#' @param cell_table a dataframe giving the number of cells per cell type
#' @import ape
#' @export
BranchCal <- function(phyla,cell_table){
  edge_list <- phyla$edge
  parent_list <- unique(edge_list[,1])
  phyla_sub <- subtrees(phyla, wait=FALSE)
  subtree_stats <- c()
  for (subtree in phyla_sub){
    tips <- subtree$tip.label
    root <- subtree$node.label[1]
    Ncells_subtree <- sum(cell_table[cell_table$Var1 %in% tips,2])
    subtree_stats <- rbind(subtree_stats, c(root,Ncells_subtree))
  }
  p_edge <- c()
  for (i in parent_list){
    childrens <- edge_list[edge_list[,1]==i,]
    ncells <- c(0,0)
    for (i in 1:2){
      children <- childrens[i,2]
      if (children <= length(phyla$tip.label)){
        ncells[i] <- cell_table[cell_table[,1]==phyla$tip.label[children],2]
      }else{
        ncells[i] <- subtree_stats[subtree_stats[,1] == children,2]
      }
    }
    prob <- ncells/sum(ncells)
    p_edge <- rbind(p_edge,cbind(childrens,prob))
  }
  return(p_edge)
}
