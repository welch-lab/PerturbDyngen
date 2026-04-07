library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ROI)
library(ROI.plugin.highs)
library(glmnet)

source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")

true_graph_test <- readRDS("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/secondary_design/example_graph_s093.rds")

# Initialize results dataframe
results_df <- data.frame(seed = integer(), precision = numeric(), recall = numeric())

for (i in 1:500) {
  set.seed(i)  # Optional if fair_ordered also uses additional randomness

  # Use the current seed for reproducibility in fair_ordered
  my_exp <- fair_ordered(222, 90, 8, seed = i)
  exp_design <- my_exp$experiment_matrix
  colnames(exp_design) <- rownames(true_graph_test)

  exp_results <- compute_edge_reachability_matrix(exp_design, true_graph_test)
  adj_DD <- infer_edges_definite(exp_design, exp_results)
  res <- compute_precision_recall_DD(adj_DD, true_graph_test)

  results_df <- rbind(
    results_df, 
    data.frame(seed = i, precision = res$precision, recall = res$recall)
  )
}

saveRDS(results_df, "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/results_DD_222_093_90_500reps.rds")