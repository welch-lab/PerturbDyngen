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


myseed = 20255
true_graph_test <- readRDS("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/secondary_design/example_graph_s093.rds")
regulators <- rownames(true_graph_test)
exp_design       <- make_random_singleton_design(57,regulators,seed = myseed)
saveRDS(exp_design,"exp_design_1st_b5.rds")
exp_results <- compute_edge_reachability_matrix(exp_design, true_graph_test)
adj_DD         <- infer_edges_definite(exp_design, exp_results)
prob_score <- compute_prob_score_matrix(exp_design, exp_results, adj_DD)



## ============ 1) 小工具：单轮评估封装 ============
safe_compute_round_metrics <- function(exp_design_all, true_graph) {
  exp_results_round <- compute_edge_reachability_matrix(exp_design_all, true_graph)
  adj_round         <- infer_edges_definite(exp_design_all, exp_results_round)
  pr                <- compute_precision_recall_DD(adj_round, true_graph)
  list(exp_results = exp_results_round, adj = adj_round, pr = pr)
}

## ============ 2) 主流程：给定 n_repeats，做二轮与三轮 ============
run_two_and_three_rounds <- function(n_repeats,
                                     prob_score_1st,
                                     adj_1st,
                                     base_exp_design,
                                     true_graph,
                                     max_size = 8,
                                     split = FALSE,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # -------- Round 2 --------
  res2 <- build_prob_based_design(
    prob_score_matrix = prob_score_1st,
    regulators        = rownames(prob_score_1st),
    adj_graph         = adj_1st,
    all_regulators    = colnames(prob_score_1st),
    max_size          = max_size,
    n_repeats         = n_repeats,
    split             = split
  )
  exp_design_2nd_all <- rbind(base_exp_design, res2$design_matrix)

  met2 <- safe_compute_round_metrics(exp_design_2nd_all, true_graph)
  pr2  <- met2$pr

  row2 <- data.frame(
    round             = "2nd",
    iterations        = n_repeats,
    n_new_experiments = nrow(res2$design_matrix),
    total_score       = res2$total_score,
    precision         = as.numeric(pr2$precision),
    recall            = as.numeric(pr2$recall),
    TP                = as.integer(pr2$TP),
    FP                = as.integer(pr2$FP),
    FN                = as.integer(pr2$FN),
    stringsAsFactors  = FALSE
  )

  # -------- Round 3（关键：重算 prob_score）--------
  prob_score_2nd <- compute_prob_score_matrix(
    exp_design  = exp_design_2nd_all,
    exp_results = met2$exp_results,
    adj_graph   = met2$adj
  )

  res3 <- build_prob_based_design(
    prob_score_matrix = prob_score_2nd,
    regulators        = rownames(prob_score_2nd),
    adj_graph         = met2$adj,
    all_regulators    = colnames(prob_score_2nd),
    max_size          = max_size,
    n_repeats         = n_repeats,
    split             = split
  )

  exp_design_3rd_all <- rbind(exp_design_2nd_all, res3$design_matrix)
  met3 <- safe_compute_round_metrics(exp_design_3rd_all, true_graph)
  pr3  <- met3$pr

  row3 <- data.frame(
    round             = "3rd",
    iterations        = n_repeats,
    n_new_experiments = nrow(res3$design_matrix),
    total_score       = res3$total_score,
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
      exp_design_2nd_add = res2$design_matrix,
      exp_design_3rd_add = res3$design_matrix
    )
  )
}

## ============ 3) 批量运行（单线程） ============
iterations_grid <- c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

out_list <- vector("list", length(iterations_grid))
for (i in seq_along(iterations_grid)) {
  iters <- iterations_grid[i]
  message(sprintf("[INFO] Running n_repeats = %d ...", iters))
  out_list[[i]] <- tryCatch(
    run_two_and_three_rounds(
      n_repeats       = iters,
      prob_score_1st  = prob_score,
      adj_1st         = adj_DD,
      base_exp_design = exp_design,
      true_graph      = true_graph_test,
      max_size        = 8,
      split           = FALSE,
      seed            = myseed
    ),
    error = function(e) {
      warning(sprintf("n_repeats=%d error: %s", iters, e$message))
      NULL
    }
  )
}

## ============ 4) 汇总与保存 ============
saveRDS(out_list, file = "res_list_b5.rds")
valid_idx  <- vapply(out_list, function(x) !is.null(x), logical(1))
out_list   <- out_list[valid_idx]
                     

metrics_df <- do.call(rbind, lapply(out_list, `[[`, "rows"))
row.names(metrics_df) <- NULL
metrics_df <- metrics_df[order(metrics_df$round, metrics_df$iterations), ]

# 保存指标
saveRDS(metrics_df, file = "prob_design_benchmark_metrics_b5.rds")
write.csv(metrics_df, file = "prob_design_benchmark_metrics_b5.csv", row.names = FALSE)

# 保存每个 n_repeats 对应的新增设计（便于复现）
artifacts <- lapply(out_list, `[[`, "artifacts")
saveRDS(artifacts, file = "prob_design_benchmark_artifacts_b5.rds")

print(metrics_df)