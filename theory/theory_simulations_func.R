
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(progress)
library(ROI)
library(ROI.plugin.highs)
library(binpackr)



construct_perturbations <- function(adj_mat, kmax = 8, fixed_size = TRUE, max_loops = 100, seed = 123) {
  set.seed(seed)
  
  num_nodes <- nrow(adj_mat)
  all_nodes <- seq_len(num_nodes)
  
  node_pairs <- combn(all_nodes, 2, simplify = FALSE)
  experiments <- list()
  
  # Helper to add unique experiment
  add_experiment <- function(new_exp) {
    new_exp <- sort(unique(new_exp))
    for (e in experiments) {
      if (identical(sort(e), new_exp)) return(FALSE)
    }
    experiments[[length(experiments) + 1]] <<- new_exp
    return(TRUE)
  }
  
  for (pair in node_pairs) {
    i <- pair[1]
    j <- pair[2]
    
    # i ∈ perturbed, j ∉ perturbed
    found <- FALSE
    for (loop in seq_len(max_loops)) {
      candidate <- if (fixed_size) {
        sample(setdiff(all_nodes, j), kmax)
      } else {
        size <- sample(1:kmax, 1)
        sample(setdiff(all_nodes, j), size)
      }
      if (i %in% candidate && !(j %in% candidate)) {
        add_experiment(candidate)
        found <- TRUE
        break
      }
    }
    if (!found) warning(sprintf("Fail: i=%d in, j=%d out", i, j))
    
    # j ∈ perturbed, i ∉ perturbed
    found <- FALSE
    for (loop in seq_len(max_loops)) {
      candidate <- if (fixed_size) {
        sample(setdiff(all_nodes, i), kmax)
      } else {
        size <- sample(1:kmax, 1)
        sample(setdiff(all_nodes, i), size)
      }
      if (j %in% candidate && !(i %in% candidate)) {
        add_experiment(candidate)
        found <- TRUE
        break
      }
    }
    if (!found) warning(sprintf("Fail: j=%d in, i=%d out", j, i))
  }
  
  # Limit experiment number < num_nodes
  if (length(experiments) >= num_nodes) {
    experiments <- experiments[seq_len(num_nodes - 1)]
  }

  # Convert to binary matrix
  exp_matrix <- matrix(0, nrow = length(experiments), ncol = num_nodes)
  for (i in seq_along(experiments)) {
    exp_matrix[i, experiments[[i]]] <- 1
  }
  rownames(exp_matrix) <- paste0("Exp", seq_len(nrow(exp_matrix)))
  colnames(exp_matrix) <- paste0("Node", all_nodes)
  
  ## -- Check ordered paired condition --
  uncovered_pairs <- list()
  for (pair in node_pairs) {
    i <- pair[1]
    j <- pair[2]
    
    # Check both directions
    i_only <- any(exp_matrix[, i] == 1 & exp_matrix[, j] == 0)
    j_only <- any(exp_matrix[, j] == 1 & exp_matrix[, i] == 0)
    
    if (!(i_only && j_only)) {
      uncovered_pairs[[length(uncovered_pairs) + 1]] <- c(i, j)
    }
  }

  return(list(
    experiment_matrix = exp_matrix,
    uncovered_pairs = uncovered_pairs,
    covered_percent = 1 - length(uncovered_pairs)/length(node_pairs)
  ))
}


generate_antichain <- function(n, m, subset_size = floor(n / 2), seed = 123) {
  set.seed(seed)
  
  full_layer <- combn(n, subset_size, simplify = FALSE)
  
  if (m > length(full_layer)) {
    stop(sprintf("Cannot generate %d antichain sets of size %d from layer %d (only %d available).",
                 m, n, subset_size, length(full_layer)))
  }
  
  sampled_sets <- sample(full_layer, m)
  return(sampled_sets)
}




construct_experiments_from_antichain <- function(index_sets, num_targets = NULL) {
  num_sets <- length(index_sets)
  all_targets <- sort(unique(unlist(index_sets)))
  if (!is.null(num_targets)) {
    all_targets <- sort(unique(c(all_targets, seq_len(num_targets))))
  }
  n <- length(all_targets)
  names(all_targets) <- paste0("J", all_targets)
  
  # Step 1: Build incidence matrix
  incidence <- matrix(0, nrow = num_sets, ncol = n)
  rownames(incidence) <- paste0("I", seq_len(num_sets))
  colnames(incidence) <- paste0("J", all_targets)
  
  for (i in seq_along(index_sets)) {
    incidence[i, paste0("J", index_sets[[i]])] <- 1
  }
  
  # Step 2: Build experiments: each Jj corresponds to a list of Iis that include it
  experiment_list <- lapply(seq_len(n), function(j) {
    which(incidence[, j] == 1)
  })
  names(experiment_list) <- paste0("J", all_targets)
  
  # Step 3: Build binary matrix: rows = experiments, columns = I's
  exp_matrix <- matrix(0, nrow = n, ncol = num_sets)
  rownames(exp_matrix) <- names(experiment_list)
  colnames(exp_matrix) <- paste0("x", seq_len(num_sets))
  
  for (j in seq_len(n)) {
    exp_matrix[j, experiment_list[[j]]] <- 1
  }

  # Step 4: Check ordered paired condition
  uncovered_pairs <- list()
  pairs <- combn(seq_len(num_sets), 2, simplify = FALSE)
  
  for (p in pairs) {
    i <- p[1]
    j <- p[2]
    i_only <- any(exp_matrix[, i] == 1 & exp_matrix[, j] == 0)
    j_only <- any(exp_matrix[, j] == 1 & exp_matrix[, i] == 0)
    
    if (!(i_only && j_only)) {
      uncovered_pairs[[length(uncovered_pairs) + 1]] <- c(i, j)
    }
  }
  
  list(
    experiment_matrix = exp_matrix,
    incidence_matrix = incidence,
    uncovered_pairs = uncovered_pairs,
    covered_percent = 1 - length(uncovered_pairs) / length(pairs)
  )
}




fair_index_sets <- function(K, l, p, max_iter = 1000) {
  # Step 1: Initial random distinct index sets
  all_combos <- combn(K, l, simplify = FALSE)
  if (length(all_combos) < p) stop("Not enough distinct index sets.")
  index_sets <- sample(all_combos, p)
  
  for (iter in seq_len(max_iter)) {
    # Count frequencies of each index from 1 to K
    flat_indices <- unlist(index_sets)
    freqs <- table(factor(flat_indices, levels = 1:K))
    
    most_freq <- as.integer(names(freqs)[which.max(freqs)])
    least_freq <- as.integer(names(freqs)[which.min(freqs)])
    
    if ((freqs[most_freq] - freqs[least_freq]) <= 1) break  # Balanced
    
    # Find a set A of size l - 1
    possible_A <- combn(setdiff(1:K, c(most_freq, least_freq)), l - 1, simplify = FALSE)
    for (A in possible_A) {
      new_M <- sort(c(most_freq, A))
      new_m <- sort(c(least_freq, A))
      if (any(sapply(index_sets, function(s) identical(s, new_M)))) {
        if (!any(sapply(index_sets, function(s) identical(s, new_m)))) {
          idx <- which(sapply(index_sets, function(s) identical(s, new_M)))[1]
          index_sets[[idx]] <- new_m
          break
        }
      }
    }
  }
  
  return(index_sets)
}



compute_l <- function(n, K) {
  # Find all l that satisfy the inequality, return the maximum one
  valid_ls <- sapply(1:K, function(l) {
    choose(K, l - 1) < n && n <= choose(K, l)
  })
  
  if (!any(valid_ls)) stop("No valid l satisfies Equation (10)")

  return(max(which(valid_ls)))  # ← max l satisfying inequality
}

fair_ordered <- function(n, K, max_nodes_per_experiment = 8, seed = 123) {
  set.seed(seed)
  
  # Step 1: Compute minimum l using Eq. 10
  l <- compute_l(n, K)
  message(sprintf("Using l = %d from Eq. (10)", l))
  
  # Step 2: Get I1,...,In from Fair(K, l, n)
  index_sets <- fair_index_sets(K, l, n)
  
  # Step 3: Build binary experiment matrix
  exp_matrix <- matrix(0, nrow = K, ncol = n)
  for (i in seq_len(n)) {
    exp_matrix[index_sets[[i]], i] <- 1
  }
  
  # Step 4: Enforce max nodes per experiment
  if (any(rowSums(exp_matrix) > max_nodes_per_experiment)) {
    stop("Generated experiments exceed max_nodes_per_experiment. Try increasing K or max.")
  }
  
  # Step 5: Check ordered pair condition
  uncovered_pairs <- list()
  all_pairs <- combn(n, 2, simplify = FALSE)
  for (pair in all_pairs) {
    i <- pair[1]
    j <- pair[2]
    i_only <- any(exp_matrix[, i] == 1 & exp_matrix[, j] == 0)
    j_only <- any(exp_matrix[, j] == 1 & exp_matrix[, i] == 0)
    if (!(i_only && j_only)) {
      uncovered_pairs[[length(uncovered_pairs) + 1]] <- c(i, j)
    }
  }
  
  list(
    experiment_matrix = exp_matrix,
    index_sets = index_sets,
    l = l,
    uncovered_pairs = uncovered_pairs,
    covered_percent = 1 - length(uncovered_pairs) / length(all_pairs)
  )
}





test_fair_ordered_min_experiments <- function(n = 1500, 
                                              K_start = 370, 
                                              step = 5, 
                                              T = 10, 
                                              max_nodes_per_experiment = 8) {
  results <- list()
  K <- K_start
  keep_going <- TRUE
  
  while (keep_going) {
    message(sprintf("Testing K = %d", K))
    successes <- 0
    designs <- list()
    
    for (t in 1:T) {
      seed <- 1000 + t
      tryCatch({
        res <- fair_ordered(n = n, 
                            K = K, 
                            max_nodes_per_experiment = max_nodes_per_experiment, 
                            seed = seed)
        successes <- successes + 1
        designs[[length(designs) + 1]] <- res$experiment_matrix
        message(sprintf("  ✓ Success (seed = %d)", seed))
      }, error = function(e) {
        message(sprintf("  ✗ Failure (seed = %d): %s", seed, e$message))
      })
    }
    
    results[[as.character(K)]] <- list(
      K = K,
      success_count = successes,
      designs = designs
    )
    
    if (successes == 0) {
      message(sprintf(" No successful designs for K = %d across %d seeds. Stopping.", K, T))
      break
    }
    
    K <- K - step
  }
  
  return(results)
}




sample_directed_graph <- function(num_tfs, num_non_tfs, sparsity = 0.1, allow_self_loop = FALSE) {
  # Total number of target nodes (columns)
  total_targets <- num_tfs + num_non_tfs
  
  # Initialize empty adjacency matrix: rows (TFs), cols (targets)
  adj_matrix <- matrix(0, nrow = num_tfs, ncol = total_targets)
  
  # All possible edges: (i,j) where i ∈ TFs and j ∈ TFs ∪ non-TFs
  possible_edges <- expand.grid(row = 1:num_tfs, col = 1:total_targets)
  
  # Exclude self-loops if not allowed (i.e., TF_i -> TF_i)
  if (!allow_self_loop) {
    possible_edges <- possible_edges[possible_edges$row != possible_edges$col, ]
  }
  
  # Determine number of edges based on sparsity
  num_edges_to_sample <- ceiling(sparsity * nrow(possible_edges))
  
  # Sample edge indices
  sampled_edges <- possible_edges[sample(nrow(possible_edges), num_edges_to_sample), ]
  
  # Fill in the adjacency matrix
  for (i in seq_len(nrow(sampled_edges))) {
    row_id <- sampled_edges$row[i]
    col_id <- sampled_edges$col[i]
    adj_matrix[row_id, col_id] <- 1
  }
  

  
  rownames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  if (num_non_tfs > 0) {

    
    #colnames(adj_matrix)[(num_tfs + 1):total_targets] <- paste0("Gene", seq_len(num_non_tfs))
    colnames(adj_matrix) <- c(paste0("TF", seq_len(num_tfs)),paste0("Gene", seq_len(num_non_tfs)))
  }
  else {
    colnames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  }

  
  return(adj_matrix)
}









compute_edge_reachability_matrix <- function(intervene_matrix, graph_matrix, error_rate = 0) {
  num_experiments <- nrow(intervene_matrix)
  num_targets <- ncol(graph_matrix)
  
  result_matrix <- matrix(0, nrow = num_experiments, ncol = num_targets)
  rownames(result_matrix) <- rownames(intervene_matrix)
  colnames(result_matrix) <- colnames(graph_matrix)
  
  for (i in seq_len(num_experiments)) {
    intervened_tfs <- which(intervene_matrix[i, ] == 1)
    
    if (length(intervened_tfs) == 0) next
    
    # Find genes reachable from any intervened TF
    reachable <- colSums(graph_matrix[intervened_tfs, , drop = FALSE]) > 0
    result_matrix[i, reachable] <- 1
    
    # Mask the intervened genes with NA
    intervened_genes <- colnames(intervene_matrix)[intervened_tfs]
    result_matrix[i, intervened_genes] <- NA
  }
  
  # Apply noise by flipping 0 <-> 1 randomly
  if (error_rate > 0) {
    non_na_mask <- !is.na(result_matrix)
    flip_candidates <- which(non_na_mask)
    num_flips <- ceiling(length(flip_candidates) * error_rate)
    
    if (num_flips > 0) {
      flip_indices <- sample(flip_candidates, num_flips)
      result_matrix[flip_indices] <- 1 - result_matrix[flip_indices]
    }
  }
  
  return(result_matrix)
}








prepare_edge_lp_matrices <- function(intervene_matrix, result_matrix, 
                                     max_out_degree = NULL, verbose = TRUE) {
  if (verbose && !requireNamespace("progress", quietly = TRUE)) {
    stop("Please install the 'progress' package with install.packages('progress')")
  }
  
  TFs <- colnames(intervene_matrix)
  Targets <- colnames(result_matrix)
  
  num_tfs <- length(TFs)
  num_targets <- length(Targets)
  
  edge_vars <- outer(TFs, Targets, paste0)
  var_names <- as.vector(edge_vars)
  var_index_map <- setNames(seq_along(var_names), var_names)
  
  constraint_matrix <- list()
  constraint_rhs <- c()
  constraint_dir <- c()
  
  total_steps <- nrow(intervene_matrix) * ncol(result_matrix)
  pb <- if (verbose) progress::progress_bar$new(
    format = "Building LP constraints [:bar] :current/:total (:percent) :elapsed",
    total = total_steps,
    width = 60,
    clear = FALSE
  ) else NULL
  
  step_counter <- 0
  
  # --- 1. Add experiment-based constraints ---
  for (e in seq_len(nrow(intervene_matrix))) {
    intervened_tfs <- which(intervene_matrix[e, ] == 1)
    tf_names <- TFs[intervened_tfs]
    
    for (t in seq_len(ncol(result_matrix))) {
      step_counter <- step_counter + 1
      if (verbose) pb$tick()
      
      outcome <- result_matrix[e, t]
      if (is.na(outcome)) next
      
      target_gene <- Targets[t]
      involved_edges <- paste0(tf_names, target_gene)
      row <- rep(0, length(var_names))
      row[var_index_map[involved_edges]] <- 1
      
      if (outcome == 1) {
        constraint_matrix[[length(constraint_matrix) + 1]] <- row
        constraint_rhs <- c(constraint_rhs, 1)
        constraint_dir <- c(constraint_dir, ">=")
      } else if (outcome == 0) {
        constraint_matrix[[length(constraint_matrix) + 1]] <- row
        constraint_rhs <- c(constraint_rhs, 0)
        constraint_dir <- c(constraint_dir, "==")
      }
    }
  }

  # --- 2. Add max out-degree constraints ---
  if (!is.null(max_out_degree)) {
    if (length(max_out_degree) == 1) {
      max_out_degree <- setNames(rep(max_out_degree, num_tfs), TFs)
    } else if (is.null(names(max_out_degree))) {
      stop("max_out_degree must be a named vector if length > 1")
    }
    
    for (tf in TFs) {
      limit <- max_out_degree[tf]
      edge_names <- paste0(tf, Targets)
      row <- rep(0, length(var_names))
      row[var_index_map[edge_names]] <- 1
      constraint_matrix[[length(constraint_matrix) + 1]] <- row
      constraint_rhs <- c(constraint_rhs, limit)
      constraint_dir <- c(constraint_dir, "<=")
    }
  }
  
  return(list(
    var_names = var_names,
    var_index_map = var_index_map,
    f.con = do.call(rbind, constraint_matrix),
    f.dir = constraint_dir,
    f.rhs = constraint_rhs
  ))
}





build_highs_model <- function(lp_data, objective = NULL) {
  n_vars <- length(lp_data$var_names)
  
  if (is.null(objective)) {
    objective <- rep(1, n_vars)
  }
  
  # Build binary MIP with ROI
  L <- OP(
    objective = objective,
    constraints = L_constraint(
      L = lp_data$f.con,
      dir = lp_data$f.dir,
      rhs = lp_data$f.rhs
    ),
    types = rep("B", n_vars),
    maximum = FALSE
  )
  
  return(L)
}


#########################################################################################


construct_structured_ordered <- function(N, k) {
  stopifnot(k <= N)
  
  g <- floor(N / k)
  r <- N %% k
  node_ids <- seq_len(N)
  
  # Step 1: Disjoint groups of size k
  disjoint_groups <- split(node_ids[1:(g * k)], rep(1:g, each = k))
  
  # Step 2: Sliding kun-style experiments
  kun_sets <- list()
  for (s in 1:(g - k + 1)) {
    group_window <- disjoint_groups[s:(s + k - 1)]
    for (t in 1:k) {
      new_exp <- sapply(group_window, function(gr) gr[t])
      kun_sets[[length(kun_sets) + 1]] <- new_exp
    }
  }
  
  # Step 3: Leftover nodes as singletons
  leftover_nodes <- node_ids[(g * k + 1):N]
  singleton_sets <- lapply(leftover_nodes, function(x) x)
  
  # Combine all experiments
  experiments <- c(disjoint_groups, kun_sets, singleton_sets)
  
  # Build binary matrix
  exp_matrix <- matrix(0, nrow = length(experiments), ncol = N)
  for (i in seq_along(experiments)) {
    exp_matrix[i, experiments[[i]]] <- 1
  }
  rownames(exp_matrix) <- paste0("Exp", seq_len(nrow(exp_matrix)))
  colnames(exp_matrix) <- paste0("Node", seq_len(N))
  
  # Check ordered pair condition
  uncovered_pairs <- list()
  all_pairs <- combn(N, 2, simplify = FALSE)
  
  for (pair in all_pairs) {
    i <- pair[1]
    j <- pair[2]
    i_only <- any(exp_matrix[, i] == 1 & exp_matrix[, j] == 0)
    j_only <- any(exp_matrix[, j] == 1 & exp_matrix[, i] == 0)
    if (!(i_only && j_only)) {
      uncovered_pairs[[length(uncovered_pairs) + 1]] <- c(i, j)
    }
  }
  
  list(
    experiment_matrix = exp_matrix,
    total_experiments = length(experiments),
    disjoint_group_count = g,
    sliding_experiments = (g - k + 1) * k,
    leftover_count = r,
    covered_percent = 1 - length(uncovered_pairs) / length(all_pairs),
    uncovered_pairs = uncovered_pairs
  )
}



############


sample_graph_exponential_degree <- function(num_tfs, num_non_tfs, rate = 2.0, allow_self_loop = FALSE, seed = 123) {
  set.seed(seed)
  total_targets <- num_tfs + num_non_tfs
  adj_matrix <- matrix(0, nrow = num_tfs, ncol = total_targets)

  for (i in 1:num_tfs) {
    k <- round(rexp(1, rate = rate))
    possible_targets <- setdiff(1:total_targets, if (!allow_self_loop) i else NULL)
    chosen <- sample(possible_targets, size = min(k, length(possible_targets)), replace = FALSE)
    adj_matrix[i, chosen] <- 1
  }

rownames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  if (num_non_tfs > 0) {
    colnames(adj_matrix)[(num_tfs + 1):total_targets] <- paste0("Gene", seq_len(num_non_tfs))
  }
  else {
    colnames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  }
  return(adj_matrix)
}




calculate_sparsity <- function(adj_matrix, ignore_na = FALSE) {
  # Flatten the matrix to a vector
  values <- as.vector(adj_matrix)
  
  if (ignore_na) {
    values <- values[!is.na(values)]
  } else {
    values[is.na(values)] <- 0  # Treat NA as zero
  }

  num_zero <- sum(values == 0)
  total_entries <- length(values)
  
  if (total_entries == 0) {
    warning("Matrix has only NA values; returning NA.")
    return(NA)
  }
  
  sparsity <- num_zero / total_entries
  return(sparsity)
}





sample_graph_geometric_degree <- function(num_tfs, num_non_tfs, prob = 0.3, allow_self_loop = FALSE) {
  total_targets <- num_tfs + num_non_tfs
  adj_matrix <- matrix(0, nrow = num_tfs, ncol = total_targets)

  for (i in 1:num_tfs) {
    k <- rgeom(1, prob = prob)
    possible_targets <- setdiff(1:total_targets, if (!allow_self_loop) i else NULL)
    chosen <- sample(possible_targets, size = min(k, length(possible_targets)), replace = FALSE)
    adj_matrix[i, chosen] <- 1
  }

rownames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  if (num_non_tfs > 0) {
    colnames(adj_matrix)[(num_tfs + 1):total_targets] <- paste0("Gene", seq_len(num_non_tfs))
  }
  else {
    colnames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  }
  return(adj_matrix)
}





sample_graph_poisson_degree <- function(num_tfs, num_non_tfs, lambda = 2.0, allow_self_loop = FALSE) {
  total_targets <- num_tfs + num_non_tfs
  adj_matrix <- matrix(0, nrow = num_tfs, ncol = total_targets)

  for (i in 1:num_tfs) {
    k <- rpois(1, lambda)
    possible_targets <- setdiff(1:total_targets, if (!allow_self_loop) i else NULL)
    chosen <- sample(possible_targets, size = min(k, length(possible_targets)), replace = FALSE)
    adj_matrix[i, chosen] <- 1
  }

rownames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  if (num_non_tfs > 0) {
    colnames(adj_matrix)[(num_tfs + 1):total_targets] <- paste0("Gene", seq_len(num_non_tfs))
  }
  else {
    colnames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  }
  return(adj_matrix)
}



sample_graph_powerlaw_degree <- function(num_tfs, num_non_tfs, gamma = 2.5, k_min = 1, allow_self_loop = FALSE, seed = 123) {
  set.seed(seed)
  total_targets <- num_tfs + num_non_tfs
  adj_matrix <- matrix(0, nrow = num_tfs, ncol = total_targets)

  rplaw <- function(n, k_min, gamma) {
    u <- runif(n)
    floor(k_min * (1 - u)^(-1 / (gamma - 1)))
  }

  degrees <- rplaw(num_tfs, k_min, gamma)
  degrees <- pmin(degrees, total_targets)

  for (i in 1:num_tfs) {
    possible_targets <- setdiff(1:total_targets, if (!allow_self_loop) i else NULL)
    chosen <- sample(possible_targets, size = min(degrees[i], length(possible_targets)), replace = FALSE)
    adj_matrix[i, chosen] <- 1
  }

rownames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  if (num_non_tfs > 0) {
    colnames(adj_matrix)[(num_tfs + 1):total_targets] <- paste0("Gene", seq_len(num_non_tfs))
  }
  else {
    colnames(adj_matrix) <- paste0("TF", seq_len(num_tfs))
  }
  return(adj_matrix)
}



calculate_precision_recall_hard <- function(true_matrix, infer_matrix) {
  infer_bin <- infer_matrix
  infer_bin[infer_bin >= 0.5] <- 1
  infer_bin[infer_bin < 0.5] <- 0

  TP <- sum(infer_bin == 1 & true_matrix == 1)
  FP <- sum(infer_bin == 1 & true_matrix == 0)
  FN <- sum(infer_bin == 0 & true_matrix == 1)

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  return(list(precision = precision, recall = recall))
}

calculate_precision_recall_zero <- function(true_matrix, infer_matrix) {
  infer_bin <- infer_matrix
  infer_bin[infer_bin > 0] <- 1
  infer_bin[infer_matrix == 0.5] <- 0

  TP <- sum(infer_bin == 1 & true_matrix == 1)
  FP <- sum(infer_bin == 1 & true_matrix == 0)
  FN <- sum(infer_bin == 0 & true_matrix == 1)

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  return(list(precision = precision, recall = recall))
}

calculate_precision_recall_soft <- function(true_matrix, infer_matrix) {
  TP <- sum(infer_matrix[true_matrix == 1])
  FP <- sum(infer_matrix[true_matrix == 0])
  FN <- sum(1 - infer_matrix[true_matrix == 1])

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  return(list(precision = precision, recall = recall))
}




calculate_T_bound <- function(E, N_TF, N_non_TF, V,k, allow_self_regulation = FALSE) {
  N_total <- N_TF * (N_TF + N_non_TF)
  if (!allow_self_regulation) {
    N_total <- N_total - N_TF
  }
  
  if (E == 0 || E >= N_total) {
    return(NA)
  }
  
  log_binom <- lgamma(N_total + 1) - lgamma(E + 1) - lgamma(N_total - E + 1)
  T_bound <- ceiling(log_binom / log(2) / (V - k))
  
  return(T_bound)
}


build_hub_groups <- function(adj, max_nodes = 8, top_k = 10, min_upstream = 2) {
  in_deg <- colSums(adj)
  hubs <- order(in_deg, decreasing = TRUE)[1:top_k]
  
  group_list <- list()
  for (h in hubs) {
    upstream <- which(adj[, h] == 1)
    if (length(upstream) < min_upstream) next  # skip if not enough upstream
    
    # If upstream > max_nodes, break into chunks
    if (length(upstream) > max_nodes) {
      chunks <- split(upstream, ceiling(seq_along(upstream) / max_nodes))
      for (ch in chunks) {
        if (length(ch) >= min_upstream) {
          group_list[[length(group_list) + 1]] <- ch
        }
      }
    } else {
      group_list[[length(group_list) + 1]] <- upstream
    }
  }
  
  return(group_list)
}




multi_pass_binpackr_exp_design <- function(priority_sets, N_tf, max_nodes = 8, max_pass = 10) {
  
  
  # Step 0: 先处理 large groups：拆分成多个 chunk
  processed_sets <- list()
  for (group in priority_sets) {
    if (length(group) > max_nodes) {
      chunks <- split(group, ceiling(seq_along(group) / max_nodes))
      processed_sets <- c(processed_sets, chunks)
    } else {
      processed_sets[[length(processed_sets) + 1]] <- group
    }
  }
  
  current_sets <- processed_sets
  
  for (pass in seq_len(max_pass)) {
    group_sizes <- sapply(current_sets, length)
    #print(group_sizes)
    
    # Step A: 不参与打包的 group（恰好 == max_nodes）：固定保留
    fixed_sets <- current_sets[group_sizes == max_nodes]
    
    # Step B: 打包 group（< max_nodes）
    small_sets <- current_sets[group_sizes < max_nodes]
    if (length(small_sets) == 0) break
    
    bin_ids <- bin_pack_ffd(sapply(small_sets, length), cap = max_nodes)
    bin_indices <- split(seq_along(small_sets), bin_ids)
    packed_sets <- lapply(bin_indices, function(idx) unlist(small_sets[idx]))
    
    new_sets <- c(fixed_sets, packed_sets)
    #print(length(new_sets ))
    
    if (length(new_sets) == length(current_sets)) {
      message(sprintf("Packing stabilized after %d passes.", pass))
      break
    }
    
    current_sets <- new_sets
  }
  
  # 构建 exp_design
  n_exp <- length(current_sets)
  exp_design <- matrix(0, nrow = n_exp, ncol = N_tf)
  for (i in seq_along(current_sets)) {
    exp_design[i, current_sets[[i]]] <- 1
  }
  
  
  return(exp_design)
}




complete_ordered_pair_design <- function(exp_design_hub, exp_design_ordered) {
  current_exp <- exp_design_hub  # [Exp × TF]
  N_tf <- ncol(current_exp)
  all_pairs <- combn(N_tf, 2, simplify = FALSE)
  
  uncovered_pairs <- list()
  for (pair in all_pairs) {
    i <- pair[1]
    j <- pair[2]
    i_only <- any(current_exp[, i] == 1 & current_exp[, j] == 0)
    j_only <- any(current_exp[, j] == 1 & current_exp[, i] == 0)
    if (!(i_only && j_only)) {
      uncovered_pairs[[length(uncovered_pairs) + 1]] <- c(i, j)
    }
  }
  
  selected_indices <- c()
  
  is_duplicate <- function(vec, mat) {
    if (nrow(mat) == 0) return(FALSE)
    apply(mat, 1, function(row) setequal(which(row == 1), which(vec == 1))) |> any()
  }
  
  for (i in seq_len(nrow(exp_design_ordered))) {
    candidate <- exp_design_ordered[i, ]
    
    # 跳过完全冗余的实验
    if (is_duplicate(candidate, current_exp)) next
    
    # 尝试加入
    current_exp <- rbind(current_exp, candidate)
    
    # 重新检查未覆盖 TF pair
    new_uncovered <- list()
    for (pair in all_pairs) {
      i <- pair[1]
      j <- pair[2]
      i_only <- any(current_exp[, i] == 1 & current_exp[, j] == 0)
      j_only <- any(current_exp[, j] == 1 & current_exp[, i] == 0)
      if (!(i_only && j_only)) {
        new_uncovered[[length(new_uncovered) + 1]] <- c(i, j)
      }
    }
    
    uncovered_pairs <- new_uncovered
    selected_indices <- c(selected_indices, i)
    
    if (length(uncovered_pairs) == 0) {
      message(sprintf("All TF pairs covered after %d extra experiments.", length(selected_indices)))
      break
    }
  }
  
  list(
    exp_design = current_exp,
    added_experiments = selected_indices,
    uncovered_pairs = uncovered_pairs,
    coverage = 1 - length(uncovered_pairs) / length(all_pairs)
  )
}




bfs_distances <- function(adj_list, start_node) {
  dist <- setNames(rep(Inf, length(adj_list)), names(adj_list))
  dist[start_node] <- 0
  queue <- c(start_node)
  
  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    
    for (neighbor in adj_list[[current]]) {
      if (dist[neighbor] == Inf) {
        dist[neighbor] <- dist[current] + 1
        queue <- c(queue, neighbor)
      }
    }
  }
  return(dist)
}




build_adj_list <- function(edges, all_nodes = NULL) {
  if (is.null(all_nodes)) {
    all_nodes <- unique(c(edges$from, edges$to))
  }
  adj_list <- setNames(vector("list", length(all_nodes)), all_nodes)
  for (i in seq_len(nrow(edges))) {
    adj_list[[edges$from[i]]] <- c(adj_list[[edges$from[i]]], edges$to[i])
  }
  return(adj_list)
}


infer_edges_definite <- function(intervene_matrix, result_matrix, verbose = TRUE) {
  if (verbose && !requireNamespace("progress", quietly = TRUE)) {
    stop("Please install the 'progress' package with install.packages('progress')")
  }
  
  p <- ncol(intervene_matrix)
  regulators <- colnames(intervene_matrix)
  targets <- colnames(result_matrix)
  
  edge_matrix <- matrix(NA, nrow = p, ncol = ncol(result_matrix),
                        dimnames = list(regulators, targets))
  
  step1_steps <- nrow(intervene_matrix) * ncol(result_matrix)
  step2_steps <- nrow(intervene_matrix) * ncol(result_matrix)
  step3_steps <- p * p
  total_steps <- step1_steps + step2_steps + step3_steps
  
  pb <- if (verbose) progress::progress_bar$new(
    format = "Inferring definite edges [:bar] :current/:total (:percent) :elapsed",
    total = total_steps,
    width = 60,
    clear = FALSE
  ) else NULL
  
  # Step 1: assign 0 if no effect observed
  for (e in seq_len(nrow(intervene_matrix))) {
    J <- which(intervene_matrix[e, ] == 1)
    for (t in seq_len(ncol(result_matrix))) {
      if (verbose) pb$tick()
      if (!is.na(result_matrix[e, t]) && result_matrix[e, t] == 0) {
        edge_matrix[J, t] <- 0
      }
    }
  }
  
  # Step 2: assign 1 if there's only one possible candidate
  for (e in seq_len(nrow(intervene_matrix))) {
    J <- which(intervene_matrix[e, ] == 1)
    for (t in seq_len(ncol(result_matrix))) {
      if (verbose) pb$tick()
      # for the current gene, if any gene in the exp is for sure a regulator, we can't make inference on more edges.
      if (any(edge_matrix[J, t] == 1, na.rm = TRUE)) next
      if (!is.na(result_matrix[e, t]) && result_matrix[e, t] == 1) {
        candidates <- J[is.na(edge_matrix[J, t])]
        if (length(candidates) == 1) {
          edge_matrix[candidates, t] <- 1
        }
      }
    }
  }
  
  # Step 3: contradiction check (skip rows with NA in result_matrix)
  for (r in seq_len(p)) {
    for (t in seq_len(p)) {
      
      if (verbose) pb$tick()
      if (r == t) next  # ignore self-loop by design
      if (!is.na(edge_matrix[r, t]) && edge_matrix[r, t] == 1) {
        valid_rows <- !is.na(result_matrix[, t])
        contradiction_found <- any(
          intervene_matrix[valid_rows, r] == 1 & result_matrix[valid_rows, t] == 0
        )
        if (contradiction_found) {
          edge_matrix[r, t] <- 9999
        }
      }
    }
    #if (verbose) pb$tick(ncol(result_matrix))
  }
  
  return(edge_matrix)
}



introduce_result_matrix_errors <- function(result_matrix, error_rate = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  logic_mask <- matrix(exp_results %in% c(0, 1), nrow = nrow(exp_results))
  
  flip_candidates <- which(logic_mask, arr.ind = TRUE)
  
  num_to_flip <- ceiling(nrow(flip_candidates) * error_rate)
  
  
  flip_indices <- flip_candidates[sample(seq_len(nrow(flip_candidates)), num_to_flip), , drop = FALSE]
  
  
  
  corrupted_result_matrix <- result_matrix
  
  for (i in seq_len(nrow(flip_indices))) {
    row <- flip_indices[i, 1]
    col <- flip_indices[i, 2]
    corrupted_result_matrix[row, col] <- 1 - corrupted_result_matrix[row, col]  # flip 0<->1
  }
  
  out_put_list <- list("err_mat" = corrupted_result_matrix, "error_indices" = flip_indices)
  
  return(out_put_list )
}



# for theoritc simulated results

compute_precision_recall_DD <- function(pred_matrix, true_matrix) {
  stopifnot(all(dim(pred_matrix) == dim(true_matrix)))
  
  # Flatten
  pred <- as.vector(pred_matrix)
  truth <- as.vector(true_matrix)
  
  # ------ Precision ------
  # Precision: Only count positions where pred is 0 or 1 (exclude NA, 9999)
  pred_valid_mask <- !is.na(pred) & pred != 9999 & !is.na(truth)
  pred_vals <- pred[pred_valid_mask]
  truth_vals <- truth[pred_valid_mask]
  
  TP <- sum(pred_vals == 1 & truth_vals == 1)
  FP <- sum(pred_vals == 1 & truth_vals == 0)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  
  # ------ Recall ------
  # Recall: for all true edges, check if they were correctly predicted (even NA is failure)
  true_edge_mask <- !is.na(truth) & truth == 1
  true_edge_pred <- pred[true_edge_mask]
  
  TP_recall <- sum(true_edge_pred == 1, na.rm = TRUE)
  FN <- sum(is.na(true_edge_pred) | true_edge_pred == 0 | true_edge_pred == 9999)
  
  recall <- if ((TP_recall + FN) > 0) TP_recall / (TP_recall + FN) else NA
  
  return(list(
    precision = precision,
    recall = recall,
    TP = TP,
    FP = FP,
    FN = FN
  ))
}



# functions used for spatialfdr in the following tasks

create_rough_graph_from_experiments <- function(exp_design, result_matrix, pval_threshold = 0.05) {
  # Step 1: Check consistency between exp_design and result_matrix
  if (!all(rownames(exp_design) == rownames(result_matrix))) {
    stop("Experiment rows in exp_design and result_matrix do not match.")
  }
  
  if (!all(colnames(exp_design) == colnames(result_matrix))) {
    stop("Gene columns in exp_design and result_matrix do not match.")
  }
  
  # Step 2: Check for inconsistent number of perturbed genes across experiments
  experiment_perturbed_genes <- rowSums(exp_design == 1)  # Sum of 1's for each experiment
  inconsistent_perturbed_genes <- table(experiment_perturbed_genes)  # Count of different numbers of perturbed genes
  
  if (length(inconsistent_perturbed_genes) > 1) {
    warning("Different number of perturbed genes across experiments. Counts: ")
    print(inconsistent_perturbed_genes)
  }
  
  # Step 3: Create edges from perturbed genes to significant genes
  # Identify significant genes (p-value < pval_threshold)
  significant_genes <- result_matrix < pval_threshold
  
  edges <- list()  # Initialize list to store edges
  
  # Iterate over each experiment
  for (i in 1:nrow(exp_design)) {
    perturbed_genes <- colnames(exp_design)[exp_design[i, ] == 1]  # Genes perturbed in this experiment
    significant_genes_in_exp <- significant_genes[i, ]  # Significant genes for this experiment
    
    # Add edges from perturbed genes to significant genes
    for (perturbed_gene in perturbed_genes) {
      for (sig_gene in names(significant_genes_in_exp[significant_genes_in_exp])) {
        if (significant_genes_in_exp[sig_gene] == TRUE && perturbed_gene != sig_gene) {
          edges <- append(edges, list(c(perturbed_gene, sig_gene)))  # Add edge (from, to)
        }
      }
    }
  }
  
  # Return the list of edges
  return(edges)
}



generate_baseline_design <- function(regulators, all_regulators, max_size = 8) {
  # ✅ Assert all regulators are valid
  if (!all(regulators %in% all_regulators)) {
    missing <- setdiff(regulators, all_regulators)
    stop("The following regulators are not in all_regulators:\n",
         paste(missing, collapse = ", "))
  }
  
  shuffled <- sample(regulators)
  n <- length(shuffled)
  
  n_full <- floor(n / max_size)
  full_chunks <- split(shuffled[1:(n_full * max_size)],
                       rep(1:n_full, each = max_size))
  
  # Proper handling when there's no leftover
  if (n_full * max_size < n) {
    leftovers <- shuffled[(n_full * max_size + 1):n]
    design_list <- c(full_chunks, as.list(leftovers))
  } else {
    design_list <- full_chunks
  }
  
  # Build binary matrix with consistent columns
  exp_mat <- matrix(0, nrow = length(design_list), ncol = length(all_regulators),
                    dimnames = list(paste0("exp", seq_along(design_list)), all_regulators))
  
  for (i in seq_along(design_list)) {
    genes <- design_list[[i]]
    if (length(genes) > 0) {
      exp_mat[i, genes] <- 1
    }
  }
  
  return(exp_mat)
}



build_score_based_design <- function(score_matrix, regulators,
                                     all_regulators,
                                     max_size = 8,
                                     n_repeats = 100) {
  best_design <- NULL
  best_score <- Inf
  
  for (trial in 1:n_repeats) {
    remaining <- sample(regulators)  # randomize input order
    design <- list()
    total_score <- 0
    
    while (length(remaining) >= max_size) {
      total <- length(remaining)
      idx <- sample(total)
      n_full <- floor(total / max_size)
      full_chunk_ids <- split(idx[1:(n_full * max_size)],
                              rep(1:n_full, each = max_size))
      
      chunks <- lapply(full_chunk_ids, function(i) remaining[i])
      
      chunk_scores <- sapply(chunks, function(set) {
        if (length(set) < 2) return(0)
        pairs <- combn(set, 2)
        sum(apply(pairs, 2, function(pair) {
          max(score_matrix[pair[1], pair[2]],
              score_matrix[pair[2], pair[1]],
              na.rm = TRUE)
        }))
      })
      
      if (any(chunk_scores == 0)) {
        best_chunks <- chunks[chunk_scores == 0]
        best_scores <- rep(0, length(best_chunks))
      } else {
        min_score <- min(chunk_scores)
        best_chunks <- chunks[chunk_scores == min_score]
        best_scores <- rep(min_score, length(best_chunks))
      }
      
      selected_idx <- sample(seq_along(best_chunks), 1)
      chosen <- best_chunks[[selected_idx]]
      chunk_score <- best_scores[[selected_idx]]
      
      design[[length(design) + 1]] <- chosen
      total_score <- total_score + chunk_score
      remaining <- setdiff(remaining, chosen)
    }
    
    if (length(remaining) > 0) {
      for (g in remaining) {
        design[[length(design) + 1]] <- g
      }
    }
    
    # Ensure consistent matrix with columns = all_regulators
    exp_mat <- matrix(0, nrow = length(design), ncol = length(all_regulators),
                      dimnames = list(paste0("exp", seq_along(design)), all_regulators))
    for (i in seq_along(design)) {
      exp_mat[i, design[[i]]] <- 1
    }
    
    if (total_score < best_score) {
      best_score <- total_score
      best_design <- exp_mat
    }
  }
  
  return(list(
    design_matrix = best_design,  # full matrix with aligned columns
    total_score = best_score
  ))
}



convert_to_design_matrix <- function(gene_sets, all_genes) {
  mat <- matrix(0, nrow = length(gene_sets), ncol = length(all_genes),
                dimnames = list(paste0("exp", seq_along(gene_sets)), all_genes))
  for (i in seq_along(gene_sets)) {
    mat[i, gene_sets[[i]]] <- 1
  }
  return(mat)
}


compute_experiment_scores <- function(exp_mat, score_matrix) {
  gene_names <- colnames(exp_mat)
  
  scores <- apply(exp_mat, 1, function(row) {
    genes <- gene_names[which(row == 1)]
    if (length(genes) < 2) return(0)
    
    pair_scores <- combn(genes, 2, function(pair) {
      max(score_matrix[pair[1], pair[2]], score_matrix[pair[2], pair[1]], na.rm = TRUE)
    })
    
    sum(pair_scores)
  })
  
  names(scores) <- rownames(exp_mat)
  return(scores)
}



compute_expected_positive_matrix <- function(exp_design, exp_results, adj) {
  genes <- rownames(adj)
  n_genes <- length(genes)
  
  score_matrix <- matrix(0, nrow = n_genes, ncol = n_genes,
                         dimnames = list(genes, genes))
  
  for (R in genes) {
    uncertain_targets <- genes[is.na(adj[R, ]) & genes != R]
    
    for (G in uncertain_targets) {
      for (R2 in setdiff(genes, c(R, G))) {
        # Condition 1: find experiments where both R and R2 were perturbed,this also helps detect internal relationships 
        # within each perturbation sets, as this penalize the genes in the same experiments in the first set
        shared_experiments <- which(exp_design[, R] == 1 & exp_design[, R2] == 1)
        # Condition 2: check if R2 is already a confirmed direct parent of G
        is_known_parent <- !is.na(adj[R2, G]) && adj[R2, G] == 1
        
        if (length(shared_experiments) > 0 || is_known_parent) {
          score_matrix[R, R2] <- score_matrix[R, R2] + 1
        }
      }
    }
  }
  
  return(score_matrix)
}

compute_prob_score_matrix <- function(exp_design, exp_results, adj_graph) {
  regulators <- rownames(adj_graph)
  targets <- colnames(adj_graph)
  
  score_matrix <- matrix(NA, nrow = length(regulators), ncol = length(targets),
                         dimnames = list(regulators, targets))
  
  for (target in targets) {
    # Find experiments where this target was tested (i.e. not perturbed)
    valid_exps <- which(!is.na(exp_results[, target]))
    
    for (reg in regulators) {
      if (!is.na(adj_graph[reg, target]) && adj_graph[reg, target] == 1) {
        score_matrix[reg, target] <- 0  # confirmed regulator
        next
      }
      
      probs <- c()
      for (exp_id in valid_exps) {
        perturbed_genes <- names(which(exp_design[exp_id, ] == 1))
        if (!(reg %in% perturbed_genes)) next  # regulator not perturbed in this exp
        
        if (exp_results[exp_id, target] != 1) next  # only consider positive results
        
        # Among perturbed genes, how many are NOT parents
        known_non_parents <- setdiff(perturbed_genes, reg)
        known_non_parents <- known_non_parents[!is.na(adj_graph[known_non_parents, target]) &
                                                 adj_graph[known_non_parents, target] == 0]
        
        total_size <- length(perturbed_genes)
        unknown_size <- total_size - length(known_non_parents)
        
        if (unknown_size <= 0) next
        
        # Denominator: total possible subsets excluding null set = 2^k - 1
        # Numerator: subsets without the given regulator = 2^(k-1) - 1
        k <- unknown_size
        prob <- (2^(k - 1) - 1) / (2^k - 1)
        
        probs <- c(probs, prob)
      }
      
      if (length(probs) == 0) {
        score_matrix[reg, target] <- 1  # no info, assume maximum uncertainty
      } else {
        score_matrix[reg, target] <- min(probs)  # most optimistic view
      }
    }
  }
  
  return(score_matrix)
}


evaluate_design_prob_score <- function(exp_design, prob_score_matrix, adj_graph) {
  regulators <- rownames(adj_graph)
  targets <- colnames(adj_graph)
  
  # Uncertain edges: NA in adj_graph
  uncertain_edges <- which(is.na(adj_graph), arr.ind = TRUE)
  total_uncertain <- nrow(uncertain_edges)
  
  edge_scores <- numeric(total_uncertain)
  
  for (i in seq_len(total_uncertain)) {
    R <- rownames(adj_graph)[uncertain_edges[i, 1]]
    G <- colnames(adj_graph)[uncertain_edges[i, 2]]
    
    if (!(R %in% colnames(exp_design))) next  # skip if regulator not in design
    
    # Experiments that perturb R
    exp_indices <- which(exp_design[, R] == 1)
    if (length(exp_indices) == 0) {
      edge_scores[i] <- 0  # no info → maximum uncertainty → never find infomrmation for these edges
      next
    }
    
    # For those experiments, take min score for this edge
    per_exp_scores <- sapply(exp_indices, function(exp_idx) {
      perturbed_genes <- colnames(exp_design)[exp_design[exp_idx, ] == 1]
      valid_genes <- intersect(perturbed_genes, rownames(prob_score_matrix))  # 安全处理
      if (length(valid_genes) == 0) return(NA)
      scores <- prob_score_matrix[valid_genes, G]
      if (all(is.na(scores))) return(NA)
      return(prod(scores, na.rm = TRUE))
    })
    
    
    # how does this edge been identified negative across all experiments.
    if (all(is.na(per_exp_scores))) {
      edge_scores[i] <- 0
    } else {
      edge_scores[i] <- max(per_exp_scores, na.rm = TRUE)
    }
  }
  
  # Final score: average over all uncertain edges
  #final_score <- mean(edge_scores, na.rm = TRUE)
  final_score <- sum(edge_scores, na.rm = TRUE)
  
  return(final_score)
}


build_prob_based_design <- function(prob_score_matrix,
                                    regulators,
                                    adj_graph,
                                    all_regulators,
                                    max_size = 8,
                                    n_repeats = 100,
                                    split = TRUE) {
  best_design <- NULL
  best_score <- -Inf  # now maximizing
  
  for (trial in 1:n_repeats) {
    shuffled <- sample(regulators)
    n <- length(shuffled)
    
    # Step 1: create full-size chunks
    n_full <- floor(n / max_size)
    full_chunks <- split(shuffled[1:(n_full * max_size)],
                         rep(1:n_full, each = max_size))
    
    # Step 2: remaining genes
    leftovers <- shuffled[(n_full * max_size + 1):n]
    
    # Split remaining into singletons or group as one
    if (length(leftovers) > 0) {
      if (split) {
        leftover_chunks <- as.list(leftovers)
      } else {
        leftover_chunks <- list(leftovers)
      }
    } else {
      leftover_chunks <- list()
    }
    
    # Step 3: combine all chunks
    design_list <- c(full_chunks, leftover_chunks)
    
    # Step 4: create binary experiment design matrix
    exp_mat <- matrix(0, nrow = length(design_list), ncol = length(all_regulators),
                      dimnames = list(paste0("exp", seq_along(design_list)), all_regulators))
    
    for (i in seq_along(design_list)) {
      exp_mat[i, design_list[[i]]] <- 1
    }
    
    # Step 5: evaluate design using probabilistic score
    total_score <- evaluate_design_prob_score(exp_mat, prob_score_matrix, adj_graph)
    
    # Step 6: keep best
    if (total_score > best_score) {
      best_score <- total_score
      best_design <- exp_mat
    }
  }
  
  return(list(
    design_matrix = best_design,
    total_score = best_score
  ))
}


weighted_bh <- function(p, w) {
  stopifnot(length(p) == length(w))
  keep <- which(!is.na(p) & !is.na(w) & w > 0)
  padj <- rep(NA_real_, length(p))
  if (length(keep) == 0) return(padj)
  
  p_use <- p[keep]
  w_use <- w[keep]
  
  # 与您现有实现一致：按 p 升序，累计权重
  o <- order(p_use)
  p_sorted <- p_use[o]
  w_sorted <- w_use[o]
  
  adj_sorted <- rev(cummin(rev(sum(w_sorted) * p_sorted / cumsum(w_sorted))))
  adj_sorted <- pmin(adj_sorted, 1)
  padj[keep[o]] <- adj_sorted
  padj
}



apply_degree_weighted_fdr_from_adj <- function(
    adj_pval,                             # matrix: regulators x targets (p-values)
    degree_mode = c("count_sig", "count_non_na"),
    pval_threshold = 0.05,                # 用于计算入度（当 degree_mode == "count_sig"）
    alpha = -1,                           # 入度→权重的幂；>0=密度越大越宽松，<0=越保守
    epsilon = 1e-8,                       # 防 0
    normalize = TRUE,                     # 归一化：sum(w)=m（m 为被检验的假设数）
    by_target = FALSE,                    # TRUE=列内做 BH；FALSE=全局一次
    test_scope = c("all", "sig_only")    # all=对所有非 NA p 做 FDR；sig_only=仅对“预显著”做 FDR
    
) {
  degree_mode <- match.arg(degree_mode)
  test_scope  <- match.arg(test_scope)
  stopifnot(is.matrix(adj_pval))
  regs <- rownames(adj_pval); tars <- colnames(adj_pval)
  
  ## 1) 计算每个 target 的候选入度
  deg_t <- sapply(seq_len(ncol(adj_pval)), function(j) {
    colp <- adj_pval[, j]
    if (degree_mode == "count_sig") {
      sum(!is.na(colp) & (colp < pval_threshold))
    } else {
      sum(!is.na(colp))
    }
  })
  names(deg_t) <- colnames(adj_pval)
  
  ## 2) target 入度 → 边级权重（同列同权；仅在该列非 NA 处赋值）
  w_mtx <- matrix(NA_real_, nrow = nrow(adj_pval), ncol = ncol(adj_pval),
                  dimnames = dimnames(adj_pval))
  for (j in seq_len(ncol(adj_pval))) {
    w_col <- (deg_t[j] + epsilon)^alpha
    w_mtx[, j] <- ifelse(!is.na(adj_pval[, j]), w_col, NA_real_)
  }
  
  ## 3) 确定“参与检验”的掩码（mask）
  # base：非 NA 的 p 都可候选
  mask <- !is.na(adj_pval)
  
  if (test_scope == "sig_only") {
    # 仅保留“预显著”的单元进入 FDR（其余单元不做校正，padj 置为 NA）
    mask <- mask & (adj_pval < pval_threshold )
  }
  # 若一个 target 列在 mask 下没有任何单元，将在后续保持 NA（不参与）
  
  # 同步把不检验单元的权重置为 NA（与 p 一致）
  w_mtx[!mask] <- NA_real_
  
  ## 4) 权重归一化仅在“将要检验”的单元上进行
  if (normalize) {
    valid_w <- w_mtx[!is.na(w_mtx) & w_mtx > 0]
    m <- sum(mask)  # 实际参与检验的假设数
    if (length(valid_w) > 0 && m > 0) {
      scale_factor <- m / sum(valid_w)
      w_mtx[!is.na(w_mtx)] <- w_mtx[!is.na(w_mtx)] * scale_factor
    }
  }
  
  ## 5) 做 weighted BH（仅在 mask==TRUE 的位置）
  padj_mtx <- matrix(NA_real_, nrow = nrow(adj_pval), ncol = ncol(adj_pval),
                     dimnames = dimnames(adj_pval))
  
  if (!by_target) {
    # 全局一次：把未检验位置的 p 与 w 置 NA 后直接跑
    p_vec <- as.vector(ifelse(mask, adj_pval, NA_real_))
    w_vec <- as.vector(w_mtx)
    padj_vec <- weighted_bh(p_vec, w_vec)
    padj_mtx[] <- matrix(padj_vec, nrow = nrow(adj_pval), ncol = ncol(adj_pval))
  } else {
    # 按列做：每列仅在 mask[, j] 为 TRUE 的位置运行
    for (j in seq_len(ncol(adj_pval))) {
      p_col <- adj_pval[, j]
      w_col <- w_mtx[, j]
      test_idx <- which(mask[, j])
      if (length(test_idx) > 0) {
        p_tmp <- rep(NA_real_, length(p_col)); p_tmp[test_idx] <- p_col[test_idx]
        w_tmp <- rep(NA_real_, length(w_col)); w_tmp[test_idx] <- w_col[test_idx]
        padj_mtx[, j] <- weighted_bh(p_tmp, w_tmp)
      } else {
        padj_mtx[, j] <- NA_real_
      }
    }
  }
  
  ## 6) 返回
  list(
    padj = padj_mtx,
    weights = w_mtx,
    target_degree = deg_t,
    tested_mask = mask,
    params = list(
      degree_mode = degree_mode,
      pval_threshold = pval_threshold,
      alpha = alpha,
      normalize = normalize,
      by_target = by_target,
      test_scope = test_scope
    )
  )
}






# Greedy 2nd-round k-mer planner (target-first, single candidate pool)
# with optional progress bar / messages
# -------------------------------------------------------------------
# Inputs:
#   adj_DD     : matrix (nTF x nTarget), values in {0,1,NA}; NA = unknown
#   exp_design : matrix (nExp x nTF),   values in {0,1}; colnames are TF names (reused verbatim)
#   result     : matrix (nExp x nTarget), values in {NA,0,1}; only 1 is positive
#   k          : integer, group size (e.g., 8)
#   progress   : one of c("bar","message","none")
#
# Returns:
#   $groups        : list of character vectors (TF names) per new experiment
#   $design_matrix : integer matrix (#groups x nTF) with 0/1; colnames identical to exp_design; no rownames
#   $score_detail  : data.frame(group, focus_target, score)
#   $total_score   : numeric
#
design_second_round_greedy_target_first <- function(adj_DD,
                                                    exp_design,
                                                    result,
                                                    k = 8,
                                                    progress = c("bar","message","none")) {
  progress <- match.arg(progress)
  
  stopifnot(is.matrix(adj_DD), is.matrix(exp_design), is.matrix(result))
  nTF     <- nrow(adj_DD)
  nTarget <- ncol(adj_DD)
  
  if (ncol(exp_design) != nTF) {
    stop("ncol(exp_design) must equal nrow(adj_DD) (TF dimension mismatch).")
  }
  if (ncol(result) != nTarget) {
    stop("ncol(result) must equal ncol(adj_DD) (target dimension mismatch).")
  }
  
  tf_names <- colnames(exp_design)  # authoritative TF names for output columns
  
  # Align adj_DD rows to tf_names if possible
  if (!is.null(rownames(adj_DD)) && !identical(rownames(adj_DD), tf_names)) {
    if (all(tf_names %in% rownames(adj_DD))) {
      adj_DD <- adj_DD[tf_names, , drop = FALSE]
    } else {
      warning("adj_DD rownames do not match exp_design colnames; proceeding by index.")
    }
  }
  
  # Masks
  U <- is.na(adj_DD)            # unknown edges per TF x target
  K <- (adj_DD == 1L)           # known edges
  results_pos <- (result == 1)  # round-1 positives
  
  # Historical positives per TF: P[i, t] = TRUE if any round-1 experiment containing TF i was positive at target t
  P <- matrix(FALSE, nTF, nTarget)
  for (i in seq_len(nTF)) {
    rows_i <- which(exp_design[, i] == 1L)
    if (length(rows_i) > 0L) {
      P[i, ] <- colSums(results_pos[rows_i, , drop = FALSE], na.rm = TRUE) > 0L
    }
  }
  C <- K | P  # conflict mask: known-edge OR historical positive
  
  # TFs that still need a second-round perturbation (have any NA edges)
  U_row_sum <- rowSums(U)
  need_idx  <- which(U_row_sum > 1L)
  if (length(need_idx) == 0L) {
    message("No TF has unknown edges; nothing to plan.")
    return(list(
      groups        = list(),
      design_matrix = matrix(0L, nrow = 0, ncol = nTF, dimnames = list(NULL, tf_names)),
      score_detail  = data.frame(group = integer(), focus_target = character(), score = double()),
      total_score   = 0
    ))
  }
  
  used <- rep(FALSE, nTF)
  groups <- list()
  group_scores <- numeric()
  group_focus  <- character()
  
  # Progress setup: we track TFs scheduled so far (out of length(need_idx))
  total_to_schedule <- length(need_idx)
  scheduled_so_far  <- 0L
  pb <- NULL
  if (progress == "bar") {
    pb <- utils::txtProgressBar(min = 0, max = total_to_schedule, initial = 0, style = 3)
    on.exit({ if (!is.null(pb)) close(pb) }, add = TRUE)
  }
  
  # Lexicographic comparison for numeric vectors
  lexicographic_less <- function(a, b) {
    stopifnot(length(a) == length(b))
    for (i in seq_along(a)) {
      if (a[i] < b[i])  return(TRUE)
      if (a[i] > b[i])  return(FALSE)
    }
    FALSE
  }
  
  # Pick best TF from pool using the merged conflict score:
  # score(c) = sum( ((U_G OR unit_at_tstar)) & C[c, ] )
  pick_best_from <- function(pool, U_G, t_star) {
    if (length(pool) == 0L) return(NA_integer_)
    conf_mask <- U_G
    conf_mask[t_star] <- TRUE  # inject the focus target into the same conflict channel
    best_c <- NA_integer_
    best_key <- NULL
    for (c in pool) {
      
      score  <- sum(conf_mask & C[c, ],na.rm = TRUE)   # primary: minimize merged conflict
      tie2   <- -sum(U[c, ])              # secondary: prefer larger unknown coverage
      tie3   <- c                         # tertiary: stable index
      key    <- c(score, tie2, tie3)
      if (is.null(best_key) || lexicographic_less(key, best_key)) {
        best_key <- key; best_c <- c
      }
    }
    best_c
  }
  
  # Main loop until all needed TFs are scheduled
  group_no <- 0L
  while (any(!used[need_idx])) {
    cand_tf  <- need_idx[!used[need_idx]]
    remaining <- length(cand_tf)
    
    # Choose focus target t* among remaining TFs: the target with the most unknowns
    unknown_per_t <- colSums(U[cand_tf, , drop = FALSE])
    if (all(unknown_per_t == 0L)) unknown_per_t <- colSums(U)  # safety
    t_star <- which.max(unknown_per_t)
    focus_target_name <- colnames(adj_DD)[t_star]
    
    group_no <- group_no + 1L
    G <- integer(0)
    U_G <- rep(FALSE, nTarget)
    pool <- cand_tf
    
    # If remaining <= k, we still greedily order all remaining; else we pick only k
    picks_needed <- if (remaining <= k) remaining else k
    for (pick in seq_len(picks_needed)) {
      c <- pick_best_from(pool, U_G, t_star)
      if (is.na(c)) break
      G <- c(G, c)
      used[c] <- TRUE
      U_G <- (U_G | U[c, ])
      pool <- setdiff(pool, c)
    }
    
    # Record group
    groups[[length(groups) + 1L]] <- tf_names[G]
    
    # Accumulate group score using the same merged conflict mask
    score_g <- 0
    U_pre <- rep(FALSE, nTarget)
    for (c in G) {
      conf_mask <- U_pre; conf_mask[t_star] <- TRUE
      score_g <- score_g + sum(conf_mask & C[c, ],na.rm = TRUE)
      U_pre <- (U_pre | U[c, ])
    }
    group_scores <- c(group_scores, score_g)
    group_focus  <- c(group_focus, focus_target_name)
    
    # Update progress
    scheduled_so_far <- scheduled_so_far + length(G)
    if (progress == "bar" && !is.null(pb)) {
      utils::setTxtProgressBar(pb, scheduled_so_far)
    } else if (progress == "message") {
      message(sprintf("Group %d: size=%d | focus_target=%s | score=%d | scheduled=%d/%d",
                      group_no, length(G), focus_target_name, score_g,
                      scheduled_so_far, total_to_schedule))
    }
    
    # If remaining <= k we are done (we just scheduled all remaining)
    if (remaining <= k) break
  }
  
  # Build round-2 design matrix: correct colnames, no rownames
  design_mat <- matrix(0L, nrow = length(groups), ncol = nTF,
                       dimnames = list(NULL, tf_names))
  for (r in seq_along(groups)) {
    cols <- match(groups[[r]], tf_names)
    design_mat[r, cols] <- 1L
  }
  
  score_df <- data.frame(
    group        = seq_along(groups),
    focus_target = group_focus,
    score        = group_scores,
    stringsAsFactors = FALSE
  )
  
  if (progress == "message") {
    message(sprintf("Done. Groups=%d, Total score=%d.", length(groups), sum(group_scores)))
  }
  
  list(
    groups        = groups,
    design_matrix = design_mat,
    score_detail  = score_df,
    total_score   = sum(group_scores)
  )
}




# Split the last (short) experiment into single-gene perturbations
# ----------------------------------------------------------------
# design_mat : integer/numeric matrix from out$design_matrix (rows = experiments, cols = TFs)
# k          : target group size (e.g., 8). If NULL, inferred as max(rowSums(design_mat)).
# Returns    : new design matrix with the last short row replaced by single-gene rows.
split_last_short_group <- function(design_mat, k = NULL) {
  stopifnot(is.matrix(design_mat))
  n <- nrow(design_mat)
  p <- ncol(design_mat)
  if (n == 0L) return(design_mat)
  
  # Infer k if not provided: use the maximum number of 1s in any row
  if (is.null(k)) {
    k <- max(rowSums(design_mat), na.rm = TRUE)
    if (!is.finite(k)) k <- 0L
  }
  
  last_idx <- n
  last_cnt <- sum(design_mat[last_idx, ] != 0L, na.rm = TRUE)
  
  # If the last row already has k (or more) perturbations, return as-is
  if (last_cnt >= k) {
    dimnames(design_mat) <- list(NULL, colnames(design_mat))
    storage.mode(design_mat) <- "integer"
    return(design_mat)
  }
  
  # If the last row has no 1s, drop it (nothing to split)
  if (last_cnt == 0L) {
    out <- design_mat[-last_idx, , drop = FALSE]
    dimnames(out) <- list(NULL, colnames(design_mat))
    storage.mode(out) <- "integer"
    return(out)
  }
  
  # Build one row per active TF in the last experiment
  ones <- which(design_mat[last_idx, ] != 0L)
  new_rows <- matrix(0L, nrow = length(ones), ncol = p)
  for (i in seq_along(ones)) new_rows[i, ones[i]] <- 1L
  colnames(new_rows) <- colnames(design_mat)
  
  # Combine: all but last row + the single-gene rows
  out <- rbind(design_mat[-last_idx, , drop = FALSE], new_rows)
  dimnames(out) <- list(NULL, colnames(design_mat))
  storage.mode(out) <- "integer"
  out
}







# Build a first-round singleton (one-hot) design
make_random_singleton_design <- function(n_exp, regulators, with_replacement = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- length(regulators)
  if (!with_replacement && n_exp > p) {
    stop("n_exp > number of regulators while with_replacement=FALSE.")
  }
  idx <- if (with_replacement) {
    sample.int(p, size = n_exp, replace = TRUE)
  } else {
    sample.int(p, size = n_exp, replace = FALSE)
  }
  M <- matrix(0L, nrow = n_exp, ncol = p,
              dimnames = list(NULL, regulators))
  M[cbind(seq_len(n_exp), idx)] <- 1L
  M
}

remove_nested_groups <- function(group_list) {
  is_subset <- function(a, b) all(a %in% b)
  keep <- rep(TRUE, length(group_list))
  
  for (i in seq_along(group_list)) {
    for (j in seq_along(group_list)) {
      if (i != j && is_subset(group_list[[i]], group_list[[j]])) {
        keep[i] <- FALSE
        break
      }
    }
  }
  return(group_list[keep])
}

generate_extra_disj_experiments <- function(existing_mat, k, n) {
  num_nodes <- ncol(existing_mat)
  new_experiments <- matrix(0, nrow = n, ncol = num_nodes)
  colnames(new_experiments) <- colnames(existing_mat)
  
  i <- 1
  while (i <= n) {
    perturbed <- sample(seq_len(num_nodes), k)
    candidate <- rep(0, num_nodes)
    candidate[perturbed] <- 1
    
    # check for duplicates
    if (!any(apply(existing_mat, 1, function(x) all(x == candidate))) &&
        !any(apply(new_experiments[seq_len(i - 1), , drop = FALSE], 1, function(x) all(x == candidate)))) {
      new_experiments[i, ] <- candidate
      i <- i + 1
    }
  }
  new_experiments
}


get_perturbed_nodes <- function(experiment_matrix, return_names = TRUE) {
  apply(experiment_matrix, 1, function(row) {
    perturbed <- which(row == 1)
    if (return_names && !is.null(colnames(experiment_matrix))) {
      return(colnames(experiment_matrix)[perturbed])
    } else {
      return(perturbed)
    }
  })
}


heck_duplicate_experiments <- function(experiment_matrix) {
  # Collapse each experiment row into a unique string (e.g., "10100")
  patterns <- apply(experiment_matrix, 1, paste0, collapse = "")
  
  # Find duplicates
  dup_idx <- duplicated(patterns)
  
  if (any(dup_idx)) {
    message("⚠️ Duplicate experiments found:")
    dup_rows <- which(dup_idx | duplicated(patterns, fromLast = TRUE))
    print(rownames(experiment_matrix)[dup_rows])
    return(experiment_matrix[dup_rows, , drop = FALSE])
  } else {
    message("✅ No duplicate experiments detected.")
    return(NULL)
  }
}



