#.libPaths("/nethome/xpan78/LinRace/Rlib")
library(devtools)
#devtools::load_all("./TedSim/")
library(parallel)
library(foreach)
library(TedSim)
library(DrImpute)
##load_all("/project/xpan78/LinRace-temp/")

library(doParallel)
registerDoParallel(cores=64)

# Simulate dataset using TedSim

#testing LinRace on simulated datasets
# ncells: 32, 128, 512
# mutation rate: 0.05, 0.1, 0.2,0.4
# 5 datasets each settings

Generate_profile_multi <- function(observed_counts, muts, states){
  ncells <- nrow(states)
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

sim_config <- (function() {
  # ncells_list <- c(32, 64, 128)
  ncells_list <- c(1024)
  mu_list <- c(0.05,0.1,0.15,0.2,0.25, 0.3, 0.35,0.4)
  run_list <- c(1:10)
  pd_list <- c(0, 1)
  tree_list <- c(0)
  Nchar_list <- c(16,64)
  pa_list <- c(0.4,0.8)

  configs = expand.grid(ncells_list, mu_list, run_list, pd_list, tree_list,Nchar_list,pa_list)
  names(configs) = c('ncells', 'mu', 'run', 'pd', 'tree','Nchar','pa')
  config_list = split(configs, seq(nrow(configs)))

  config_list
})()

foreach(conf = sim_config) %dopar% {
  tree_ct <- readRDS('/project/hli691/linrace/tree_ct.rds')

  #ncells <- conf$ncells
  ncells <- 1024
  mu <- conf$mu
  run <- conf$run
  p_d <- conf$pd
  tree <- conf$tree
  N_char <- conf$Nchar
  p_a <- conf$pa


  if (tree == 0) {
    phyla <- read.tree(text='((t1:4, t2:4, t3:4, t4:4, t5:4, t6:4):2);')
  } else {
    phyla <- read.tree(text='((t1:4, (t2:2, t3:2):2):2);')
  }
  #phyla <- read.tree(text='((t1:2, t2:2):1, (t3:2, t4:2):1):2;')
  N_nodes <- 2*ncells-1
  ngenes <- 500
  max_walk <- 5
  #p_a <- 0.8
  n_cif <- 30
  n_diff <- 20
  cif_step <- 0.5
  #N_char <- 16
  set.seed(run)


  returnlist <- SIFGenerate(phyla,n_diff,step = cif_step)
  cifs <- SimulateCIFs(ncells,phyla,p_a = p_a,n_CIF = n_cif,n_diff = n_diff,step = cif_step,p_d = p_d, Sigma = 0.5,mu = mu, N_char = N_char, max_walk = max_walk, SIF_res = returnlist, unif_on = FALSE)

  #We only need the leaf cells for experiments
  cif_leaves <- lapply(c(1:3),function(parami){
    cif_leaves_all <- cifs[[1]][[parami]][c(1:ncells),]
    return(cif_leaves_all)
  })
  cif_res <- list(cif_leaves,cifs[[2]])
  states <- cifs[[2]]
  states <- states[1:N_nodes,]
  states_leaves <- states[1:ncells,]
  muts <- cifs[[7]]
  rownames(muts) <- paste("cell",states[1:(N_nodes-1),4],sep = "_")
  muts_leaves <- muts[1:ncells,]

  tree_ct <- cifs[[4]]
  tip_label <- c()
  for (node in tree_ct$tip.label){
    tip <- paste('cell_',substr(node,2,nchar(node)),sep = '')
    tip_label <- c(tip_label,tip)
  }
  tree_ct$tip.label <- tip_label

  #simulate true counts
  true_counts_res <- CIF2Truecounts(ngenes = 500,ncif = n_cif,ge_prob = 0.3,ncells = ncells, cif_res = cif_res)

  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.2, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)
  profile_res <- Generate_profile_multi(observed_counts,muts,states_leaves)
  profile_out <- profile_res[[1]]
  top_genes <- profile_res[[2]]

  #gene_expression_dir <- sprintf("./sim/tree%d/expr_%d_pa_%d_mu_%d_pd_%d_Nchar_$d_run_%d.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  #cell_meta_dir <- sprintf("./sim/tree%d/meta_%d_pa_%d_mu_%d_pd_%d_Nchar_$d_run_%d.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  #tree_gt_dir <- sprintf("./sim/tree%d/tree_%d_pa_$d_mu_%g_pd_%d_Nchar_$d_run_%d.tree", tree,p_a, ncells, mu, p_d,N_char, run)
  #mut_leaves_dir <- sprintf("./sim/tree%d/mut_%d_pa_$d_mu_%g_pd_%d_Nchar_$d_run_%d.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  #cif_dir <- sprintf("./sim/tree%d/cif_%d_pa_$d_mu_%g_pd_%d_Nchar_$d_run_%d.Rds", tree,p_a, ncells, mu, p_d,N_char, run)
  #combined_profile_dir <- sprintf("./sim/tree%d/profile_%d_pa_$d_mu_%g_pd_%d_Nchar_$d_run_%d.txt", tree,p_a, ncells, mu, p_d,N_char, run)
  #top_genes_dir <- sprintf("./sim/tree%d/top_genes_%d_pa_$d_mu_%g_pd_%d_Nchar_$d_run_%d.txt", tree,p_a, ncells, mu, p_d,N_char, run)

  order <- sample(nrow(muts_leaves))
  states_leaves <- states_leaves[order,]
  muts_leaves <- muts_leaves[order,]
  counts <- true_counts_res[[1]][order,]
  profile_out <- profile_out[order,]

  gene_expression_dir <- sprintf("/project/xpan78/sim_rand/tree%d/expr_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  cell_meta_dir <- sprintf("/project/xpan78/sim_rand/tree%d/meta_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  tree_gt_dir <- sprintf("/project/xpan78/sim_rand/tree%d/tree_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.tree", tree,p_a, ncells, mu, p_d,N_char, run)
  mut_leaves_dir <- sprintf("/project/xpan78/sim_rand/tree%d/mut_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.csv", tree,p_a, ncells, mu, p_d,N_char, run)
  cif_dir <- sprintf("/project/xpan78/sim_rand/tree%d/cif_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.Rds", tree,p_a, ncells, mu, p_d,N_char, run)
  combined_profile_dir <- sprintf("/project/xpan78/sim_rand/tree%d/profile_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.txt", tree,p_a, ncells, mu, p_d,N_char, run)
  top_genes_dir <- sprintf("/project/xpan78/sim_rand/tree%d/top_genes_%g_pa_%g_mu_%g_pd_%g_Nchar_%g_run_%g.txt", tree,p_a, ncells, mu, p_d,N_char, run)

  print("writing to directory...")
  write.csv(true_counts_res[[1]], gene_expression_dir, row.names = FALSE)
  write.csv(states_leaves, cell_meta_dir)
  write.csv(muts_leaves, mut_leaves_dir)
  write.tree(cifs[[4]], tree_gt_dir)
  saveRDS(cifs, cif_dir)
  write.table(profile_out,file = combined_profile_dir,row.names = FALSE,sep = "\t", quote = FALSE)
  write.table(subset(top_genes,select = 1),top_genes_dir,row.names = FALSE,sep = "\t", quote = FALSE)

}
