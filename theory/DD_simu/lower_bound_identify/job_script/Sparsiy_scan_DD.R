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


sparsities <- c(0.2, 0.15, 0.10, 0.05, 0.01)
num_exp <- 65
my_exp <- fair_ordered(250, num_exp, 8)
results_exp2_sparsity_based <- list()

for (s in sparsities) {
  tmp_list <- replicate(5, {
    true_graph_test <- sample_directed_graph(250, 0, sparsity = s)
    exp_design <- my_exp$experiment_matrix
    colnames(exp_design) <- rownames(true_graph_test)
    exp_results <- compute_edge_reachability_matrix(exp_design, true_graph_test)
    adj_DD <- infer_edges_definite(exp_design, exp_results)

    metrics <- compute_precision_recall_DD(adj_DD, true_graph_test)
    F1 <- if (!is.na(metrics$precision) && !is.na(metrics$recall) &&
                (metrics$precision + metrics$recall) > 0) {
      2 * metrics$precision * metrics$recall / (metrics$precision + metrics$recall)
    } else NA

    list(metrics = metrics, F1 = F1)
  }, simplify = FALSE)

  avg_precision <- mean(sapply(tmp_list, function(x) x$metrics$precision), na.rm = TRUE)
  avg_recall <- mean(sapply(tmp_list, function(x) x$metrics$recall), na.rm = TRUE)
  avg_F1 <- mean(sapply(tmp_list, function(x) x$F1), na.rm = TRUE)

  results_exp2_sparsity_based[[as.character(s)]] <- data.frame(
    sparsity = s,
    metric = c("precision", "recall", "F1"),
    value = c(avg_precision, avg_recall, avg_F1)
  )
}

df_exp2_sparsity <- bind_rows(results_exp2_sparsity_based)
out_csv <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/lower_bound_identify/precision_recall_F1_vs_sparsity_N250.csv"
write.csv(df_exp2_sparsity, out_csv, row.names = FALSE)

out_png <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/DD_simu/lower_bound_identify/precision_recall_F1_vs_sparsity_N250.png"
png(out_png, width = 12, height = 8, units = "in", res = 300)
ggplot(df_exp2_sparsity, aes(x = sparsity, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Precision, Recall & F1 vs. Sparsity (DD alg, Fixed #Experiments)",
       y = "Score") +
  theme_minimal()
dev.off()
