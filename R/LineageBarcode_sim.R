#' Generate mutated states from a real dataset
#' @param N_ms the number of mutated states required
#' @param cm character matrix of the real dataset
Mutated_state_dist <- function(N_ms,cm){
  cm_table <- table(unlist(cm))
  cm_table <- sort(cm_table,decreasing = TRUE)
  cm_table <- cm_table[dimnames(cm_table)[[1]]!='0']
  cm_table <- cm_table[dimnames(cm_table)[[1]]!='-']

  if (N_ms < length(dimnames(cm_table)[[1]])){
    cm_table <- cm_table[1:N_ms]
  }
  cm_prop <- prop.table(cm_table)
  states_dist <- cm_prop
  return(states_dist)
}

#' Generate CRISPR induced mutations
#' @param barcode input barcode (from the parent cell)
#' @param mu mutation rate per target, default is 1
#' @param N_ms the number of mutated states required
#' @param mutation_dist input distribution for the mutated states
#' @param p_d whether or not to simulate dropout effects, 0 or 1
#' @param unif_on distribution mode of mutated states. TRUE: uniform distribution; FALSE: real fitted distribution
generate_mutation <- function(barcode, mu=0.1 , N_ms = 100 ,mutation_dist = NULL, p_d = 0, unif_on = FALSE){
  if (unif_on){
    states <- as.list(seq (1,N_ms,1))
    prob_dist <- NULL
  }else{
    states <- dimnames(mutation_dist)[[1]]
    prob_dist <- as.vector(mutation_dist)
  }

  m <- length(barcode)
  child_barcode <- barcode
  mu_loc = runif(m) < mu
  mutation_cites = (child_barcode == 0) &  mu_loc
  n_mut = sum(mutation_cites)
  if (n_mut != 0) {
    child_barcode[mutation_cites] = as.integer(sample(states, n_mut, replace = T, prob = prob_dist))
    if ((n_mut >=2)&(p_d == 1)){
      child_barcode <- generate_dropout(child_barcode,mutation_cites)
    }
  }

  return (child_barcode)
}

#' Generate excision dropout
#' @param barcode input barcode with mutations
#' @param mutation_cites the mutated target sites
generate_dropout <- function(barcode,mutation_cites){
  dropout_between = sample(which(mutation_cites), 2 )
  barcode[dropout_between[1]:dropout_between[2]] <- Inf
  return(barcode)
}


#' Simulate Cell Identity Factor Matrix With Lineage Barcodes
#' @param observed_counts Observed counts of gene expressions
#' @param muts Unprocessed mutation barcodes
#' @export
CaptureDrop <- function(observed_counts, muts) {
  #find the top 10 highly expressed genes in the observed counts
  #observed_counts <- t(observed_counts)
  n_char <- dim(muts)[2]
  observed_counts <-  log(observed_counts+1)
  gene_means <- rowMeans(observed_counts)
  gene_order <- order(gene_means,decreasing = TRUE)[1:10]
  reg <- sample(1:10,1)
  gene_reg <- observed_counts[,gene_order[reg]]
  drop_reg <- which(gene_reg == 0)
  N_drop <- length(drop_reg)
  for (i in drop_reg){
    muts[i,] <- rep('-',n_char)
  }
  print(sprintf("Barcodes of %d cells removed based on gene %d.", N_drop, gene_order[reg] ))
  return(muts)
}
