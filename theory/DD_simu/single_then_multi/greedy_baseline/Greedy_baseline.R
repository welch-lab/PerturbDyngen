library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
# library(lpSolve)
library(progress)
library(ROI)
library(ROI.plugin.highs)
library(glmnet)

source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")

# ---------- Inputs ----------
true_graph_test <- readRDS("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/secondary_design/example_graph_s093.rds")
regulators <- rownames(true_graph_test)
n_first_round <- 57L
k <- 8L

# ---------- Helpers ----------
safe_compute_round_metrics <- function(exp_design_all, true_graph) {
  exp_results_round <- compute_edge_reachability_matrix(exp_design_all, true_graph)
  adj_round         <- infer_edges_definite(exp_design_all, exp_results_round)
  pr                <- compute_precision_recall_DD(adj_round, true_graph)
  list(exp_results = exp_results_round, adj = adj_round, pr = pr)
}

# Greedy-only 2nd & 3rd rounds; iterations is fixed to -1; total_score from out_2nd/out_3rd
run_two_and_three_rounds <- function(base_exp_design,
                                     adj_1st,
                                     true_graph,
                                     max_size = 8,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Round 2 (greedy target-first)
  exp_results_1st <- compute_edge_reachability_matrix(base_exp_design, true_graph)
  out_2nd <- design_second_round_greedy_target_first(
    adj_DD      = adj_1st,
    exp_design  = base_exp_design,
    result = exp_results_1st,
    k           = max_size
  )
  exp_design_2nd =  split_last_short_group(out_2nd$design_matrix, k = 8)
  exp_design_2nd_all <- rbind(base_exp_design,exp_design_2nd)

  met2 <- safe_compute_round_metrics(exp_design_2nd_all, true_graph)
  pr2  <- met2$pr

  row2 <- data.frame(
    round             = "2nd",
    iterations        = -1L,
    n_new_experiments = nrow(exp_design_2nd),
    total_score       = out_2nd$total_score,
    precision         = as.numeric(pr2$precision),
    recall            = as.numeric(pr2$recall),
    TP                = as.integer(pr2$TP),
    FP                = as.integer(pr2$FP),
    FN                = as.integer(pr2$FN),
    stringsAsFactors  = FALSE
  )

  # Round 3 (greedy target-first on updated knowledge)
  exp_results_2nd <- met2$exp_results
  adj_2nd         <- met2$adj

  out_3rd <- design_second_round_greedy_target_first(
    adj_DD      = adj_2nd,
    exp_design  = exp_design_2nd_all,
    result = exp_results_2nd,
    k           = max_size
  )
  exp_design_3rd =  split_last_short_group(out_3rd$design_matrix, k = 8)
  exp_design_3rd_all <- rbind(exp_design_2nd_all, exp_design_3rd )

  met3 <- safe_compute_round_metrics(exp_design_3rd_all, true_graph)
  pr3  <- met3$pr

  row3 <- data.frame(
    round             = "3rd",
    iterations        = -1L,
    n_new_experiments = nrow(exp_design_3rd),
    total_score       = out_3rd$total_score,
    precision         = as.numeric(pr3$precision),
    recall            = as.numeric(pr3$recall),
    TP                = as.integer(pr3$TP),
    FP                = as.integer(pr3$FP),
    FN                = as.integer(pr3$FN),
    stringsAsFactors  = FALSE
  )

  list(
    rows = rbind(row2, row3),
    artifacts = list(
      exp_design_2nd_add = exp_design_2nd,
      exp_design_3rd_add = exp_design_3rd
    )
  )
}

# ---------- Main: seeds 202501:202505 ----------
seed_list <- 20251:20255
out_list  <- vector("list", length(seed_list))

for (i in seq_along(seed_list)) {
  s <- seed_list[i]
  message(sprintf("[INFO] Running greedy pipeline with seed = %d ...", s))

  # 1st-round design: random singleton without replacement (57 unique targets)
  exp_design <- make_random_singleton_design(n_first_round, regulators, seed = s)
  saveRDS(exp_design, sprintf("greedy_exp_design_1st_seed%d.rds", s))

  # 1st-round results and DD graph (for round-2 greedy input)
  exp_results <- compute_edge_reachability_matrix(exp_design, true_graph_test)
  adj_DD      <- infer_edges_definite(exp_design, exp_results)

  # Run greedy 2nd & 3rd rounds
  out_list[[i]] <- tryCatch(
    run_two_and_three_rounds(
      base_exp_design = exp_design,
      adj_1st         = adj_DD,
      true_graph      = true_graph_test,
      max_size        = k,
      seed            = s
    ),
    error = function(e) {
      warning(sprintf("seed=%d error: %s", s, e$message))
      NULL
    }
  )
}

# ---------- Collect & save ----------
saveRDS(out_list, file = "greedy_res_list.rds")
valid_idx <- vapply(out_list, function(x) !is.null(x), logical(1))
out_list  <- out_list[valid_idx]

metrics_df <- do.call(rbind, lapply(out_list, `[[`, "rows"))
row.names(metrics_df) <- NULL
metrics_df <- metrics_df[order(metrics_df$round, metrics_df$iterations), ]

saveRDS(metrics_df, file = "greedy_benchmark_metrics.rds")
write.csv(metrics_df, file = "greedy_benchmark_metrics.csv", row.names = FALSE)

artifacts <- lapply(out_list, `[[`, "artifacts")
saveRDS(artifacts, file = "greedy_benchmark_artifacts.rds")

print(metrics_df)
