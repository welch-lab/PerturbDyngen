library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
#library(lpSolve)
library(progress)
library(ROI)
library(ROI.plugin.highs)
library(glmnet)
source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")




my_exp <- fair_ordered(1500,375,8)
#true_graph_test = sample_directed_graph(1500,100000,sparsity = 0.07)
true_graph_test = sample_directed_graph(1500,100000,sparsity = 0.05)
print(sum(true_graph_test)/(dim(true_graph_test)[1] * dim(true_graph_test)[2]))
exp_design = my_exp$experiment_matrix
colnames(exp_design) = rownames(true_graph_test)
exp_results = compute_edge_reachability_matrix(exp_design,true_graph_test)
adj_DD = infer_edges_definite(exp_design, exp_results)
res = compute_precision_recall_DD(adj_DD,true_graph_test )
saveRDS(adj_DD,"/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/adj_DD_005.rds")

saveRDS(res,"/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/adj_DD_res_005.rds")