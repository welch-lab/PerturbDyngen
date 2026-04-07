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



true_graph_test <- sample_directed_graph(250, 0, sparsity = 0.01)

experiment_counts <- seq(65, 250, by = 20)
results_exp1 <- list()

for (num_exp in experiment_counts) {
  my_exp <- fair_ordered(250, num_exp, 8)
  exp_design <- my_exp$experiment_matrix
  colnames(exp_design) <- rownames(true_graph_test)

  exp_results <- compute_edge_reachability_matrix(exp_design, true_graph_test)
  adj_DD <- infer_edges_definite(exp_design, exp_results)

  metrics <- compute_precision_recall_DD(adj_DD, true_graph_test)
  F1 <- if (!is.na(metrics$precision) && !is.na(metrics$recall) &&
              (metrics$precision + metrics$recall) > 0) {
    2 * metrics$precision * metrics$recall / (metrics$precision + metrics$recall)
  } else NA

  results_exp1[[as.character(num_exp)]] <- data.frame(
    experiments = num_exp,
    metric = c("precision", "recall", "F1"),
    value = c(metrics$precision, metrics$recall, F1)
  )
}

df_exp1 <- bind_rows(results_exp1)
out_csv <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/lower_bound_identify/precision_recall_F1_vs_experiments_N250_sp99.csv"
write.csv(df_exp1, out_csv, row.names = FALSE)

# 理论下界
T_bound <- calculate_T_bound(sum(true_graph_test == 1), 250, 0, 250, 8)

out_png <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/lower_bound_identify/precision_recall_F1_vs_experiments_N250_sp99.png"
png(out_png, width = 12, height = 8, units = "in", res = 300)
ggplot(df_exp1, aes(x = experiments, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_point() +
  geom_vline(xintercept = T_bound, size = 2) +
  annotate("text", x = T_bound, y = max(df_exp1$value, na.rm = TRUE) * 0.2,
           label = paste("Lower Bound =", T_bound, 
                         ", Sparsity =", calculate_sparsity(true_graph_test)),
           vjust = -0.5, angle = 90, hjust = -0.1) +
  labs(title = "Precision, Recall & F1 vs. #Experiments (DD alg, Fixed Graph)",
       y = "Score") +
  theme_minimal()
dev.off()