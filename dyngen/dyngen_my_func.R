
compute_precision_recall <- function(adj_truth, adj_pvals, cutoff = 0.05, fill_na = TRUE, binary = TRUE) {
  # Ensure both matrices have the same nodes and order
  common_nodes <- intersect(rownames(adj_truth), rownames(adj_pvals))
  
  # Reorder both matrices to align by common nodes
  adj_truth <- adj_truth[common_nodes, common_nodes]
  adj_pvals <- adj_pvals[common_nodes, common_nodes]
  
  # Compute the binary ground truth adjacency (1 for connection, 0 otherwise)
  if (binary){
  ground_truth_bin <- sign(abs(adj_truth))}
  
  if(fill_na){
    adj_pvals[is.na(adj_pvals)] <- 1
  }
  # Compute the binary prediction based on the cutoff for p-values
  predicted_bin <- ifelse(adj_pvals < cutoff, 1, 0)
  
  # Flatten matrices for element-wise comparison
  gt_vector <- as.vector(ground_truth_bin)
  pred_vector <- as.vector(predicted_bin)
  
  # Compute true positives, false positives, and false negatives
  true_positives <- sum(gt_vector == 1 & pred_vector == 1)
  false_positives <- sum(gt_vector == 0 & pred_vector == 1)
  false_negatives <- sum(gt_vector == 1 & pred_vector == 0)
  
  # Compute precision and recall, handling division by zero cases
  precision <- ifelse((true_positives + false_positives) > 0, 
                      true_positives / (true_positives + false_positives), 
                      NA)
  recall <- ifelse((true_positives + false_negatives) > 0, 
                   true_positives / (true_positives + false_negatives), 
                   NA)
  
  # Return results as a list
  return(list(precision = precision, recall = recall, 
              true_positives = true_positives, 
              false_positives = false_positives, 
              false_negatives = false_negatives,
              adj_truth = adj_truth,
              adj_pvals = adj_pvals,
              predicted_bin = predicted_bin,
              ground_truth_bin = ground_truth_bin))
}

# not consider the perturbed cells number should be similar, this take wt and perturbed together
create_equal_count_bins <- function(sim_times, num_bins) {
  sorted_times <- sort(sim_times)
  n <- length(sorted_times)
  
  # Define breakpoints based on quantiles
  bin_edges <- quantile(sorted_times, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
  bin_edges <- unique(bin_edges)# Ensure unique edges
  # Assign each time value to a bin
  bin_labels <- cut(sim_times, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
  
  # Return a data frame
  return(data.frame(SimTime = sim_times, Bin = bin_labels))
}

create_equal_count_bins_perturb <- function(df, num_bins) {
  # Filter for perturbed data
  perturbed_df <- df[df[["Type"]] == "perturbed", ]
  sim_times_perturbed <- perturbed_df$SimTime
  
  # Compute bin edges from perturbed only
  bin_edges <- quantile(sim_times_perturbed, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
  bin_edges <- unique(bin_edges)  # Ensure unique bin edges
  
  # Assign bin labels to all rows based on these edges
  df$Bin <- cut(df$SimTime, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
  
  # Count number of observational points per bin
  obs_counts <- table(df$Bin[df$Type == "observational"])
  
  # Create a vector to match observational count per bin
  df$ObsCountPerBin <- obs_counts[as.character(df$Bin)]
  df$ObsCountPerBin[is.na(df$ObsCountPerBin)] <- 0  # Set NA (bins with no obs data) to 0
  
  return(df)
}

plot_simulation_expression_models <- function(models, simulation_i = 1:4, what = c("mol_premrna", 
    "mol_mrna", "mol_protein"), facet = c("simulation", "module_group", 
    "module_id", "none"), label_nonzero = FALSE) 
{
    is_tf <- mol <- val <- molecule <- value <- module_id <- type <- sim_time <- color <- module_group <- x <- y <- condition <- NULL
    assertthat::assert_that(what %all_in% c("mol_premrna", "mol_mrna", "mol_protein"))
    facet <- match.arg(facet)

    # Combine WT and experimental models into one dataframe
    df_list <- lapply(names(models), function(model_name) {
        model <- models[[model_name]]
        molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, !!!what) %>% pull(val)
        
        df <- bind_cols(model$simulations$meta, as.data.frame(as.matrix(model$simulations$counts)[, molecules])) %>%
            filter(simulation_i %in% !!simulation_i) %>% 
            gather(molecule, value, one_of(molecules)) %>%
            left_join(model$feature_info %>% 
                select(!!!what, module_id) %>% gather(type, molecule, !!!what) %>% 


         mutate(type = factor(type, levels = what)), by = "molecule") %>%
            group_by(module_id, sim_time, simulation_i, type) %>%
            summarise(value = mean(value), .groups = "drop") %>%
            ungroup() %>%
            mutate(module_group = gsub("[0-9]*$", "", module_id),
                   condition = model_name)  # Add model name as a condition identifier

        return(df)
    })

    df <- bind_rows(df_list)  # Combine all models into a single dataframe

    # Plot both WT and experimental conditions for comparison
    g <- ggplot(df, aes(sim_time, value, colour = condition, linetype = type)) +
        geom_step(aes(size = type)) +
        scale_size_manual(values = c(mol_premrna = 0.3, mol_mrna = 0.8, mol_protein = 0.3)) +
        scale_colour_manual(values = setNames(rainbow(length(unique(df$condition))), unique(df$condition))) + 
        theme_bw() +
        labs(colour = "Condition", linetype = "Molecule type")

    if (label_nonzero) {
        pts <- seq(0, max(df$sim_time), by = 5)
        df_labels <- df %>%
            group_by(module_id, type, module_group, condition) %>%
            do({
                df2 <- .
                approx(x = df2$sim_time, y = df2$value, xout = pts) %>%
                    as_tibble() %>%
                    rename(sim_time = x, value = y)
            }) %>%
            ungroup() %>%
            filter(value > 0)
        g <- g + geom_point(data = df_labels) + 
            ggrepel::geom_text_repel(aes(label = paste0(module_id, "_", type)), df_labels)
    }

    # Apply faceting based on the chosen parameter
    if (facet == "simulation") {
        g <- g + facet_wrap(~simulation_i, ncol = 1)
    } else if (facet == "module_group") {
        g <- g + facet_wrap(~module_group, ncol = 1)
    } else if (facet == "module_id") {
        g <- g + facet_wrap(~module_id, ncol = 1)
    }

    return(g)
}

generate_cells_debug = 
function (model) 
{
    time <- NULL
    if (model$verbose) 
        cat("Precompiling reactions for simulations\n")
    model <- dyngen:::.add_timing(model, "6_simulations", "precompile reactions for simulations")
    reactions <- dyngen:::.generate_cells_precompile_reactions(model)
    if (model$verbose) 
        cat("Running ", nrow(model$simulation_params$experiment_params), 
            " simulations\n", sep = "")
    model <- dyngen:::.add_timing(model, "6_simulations", "running simulations")
    simulations <- pbapply::pblapply(X = seq_len(nrow(model$simulation_params$experiment_params)), 
        cl = model$num_cores, FUN = dyngen:::.generate_cells_simulate_cell, 
        model = model, reactions = reactions)
    model <- dyngen:::.add_timing(model, "6_simulations", "generate output")
    model$simulations <- lst(meta = map_df(simulations, "meta"), 
        counts = do.call(rbind, map(simulations, "counts")), 
        cellwise_grn = do.call(rbind, map(simulations, "cellwise_grn")), 
        reaction_firings = do.call(rbind, map(simulations, "reaction_firings")), 
        reaction_propensities = do.call(rbind, map(simulations, 
            "reaction_propensities")), rna_velocity = do.call(rbind, 
            map(simulations, "rna_velocity")), kd_multiplier = do.call(rbind, 
            map(simulations, "kd_multiplier")), perturbed_parameters = do.call(rbind, 
            map(simulations, "perturbed_parameters")))
    if (model$verbose) 
        cat("Mapping simulations to gold standard\n", sep = "")
    model <- dyngen:::.add_timing(model, "6_simulations", "map simulations to gold standard")
    if (!is.null(model[["gold_standard"]])) {
        model$simulations$meta <- dyngen:::.generate_cells_predict_state(model)
        print(model$simulations$meta)
    }
    else {
        model$simulations$meta <- model$simulations$meta %>% 
            rename(sim_time = time)
    }
    model <- dyngen:::.add_timing(model, "6_simulations", "perform dimred")
    if (model$simulation_params$compute_dimred) {
        if (model$verbose) 
            cat("Performing dimred\n", sep = "")
        model <- model %>% dyngen:::calculate_dimred()
    }
    dyngen:::.add_timing(model, "6_simulations", "end")
}

generate_experiment_sample_cells_all = function(model){
    network <- model$gold_standard$network
    end_states <- setdiff(unique(network$to), unique(network$from))
    sim_meta <- model$simulations$meta %>% mutate(orig_ix = row_number()) %>% 
        filter(.data$sim_time >= 0)
    params <- model$experiment_params
    params$fun(network = network, sim_meta = sim_meta, params = model$experiment_params, 
        num_cells = model$numbers$num_cells)
}

generate_experiment_debug = 
function (model) 
{
    if (model$verbose) 
        cat("Simulating experiment\n")
    model <- dyngen:::.add_timing(model, "7_experiment", "sample cells")
    sample_df <- generate_experiment_sample_cells_all(model)
    if (!is.data.frame(sample_df)) {
        sample_df <- tibble(step_ix = sample_df)
    }
    cell_info <- bind_cols(model$simulations$meta[sample_df$step_ix, 
        , drop = FALSE], sample_df) %>% sample_n(n(), replace = FALSE) %>% 
        mutate(cell_id = paste0("cell", row_number())) %>% select(.data$cell_id, 
        .data$step_ix, .data$simulation_i, .data$sim_time, .data$from, 
        .data$to, .data$time, everything())
    step_ixs <- cell_info$step_ix
    tsim_counts <- model$simulations$counts[step_ixs, , drop = FALSE]
    rownames(tsim_counts) <- cell_info$cell_id
    mol_info <- model$feature_info %>% gather("mol", "val", .data$mol_premrna, 
        .data$mol_mrna, .data$mol_protein) %>% select(.data$feature_id, 
        .data$mol, .data$val)
    model <- dyngen:::.add_timing(model, "7_experiment", "fetch realcount")
    realcount <- dyngen:::.generate_experiment_fetch_realcount(model)
    model <- dyngen:::.add_timing(model, "7_experiment", "simulate library size variation")
    mrna_ids <- mol_info %>% filter(.data$mol != "mol_protein") %>% 
        pull(.data$val)
    tsim_counts_mrna <- tsim_counts[, mrna_ids, drop = FALSE]
    count_mrna_simulation <- dyngen:::.simulate_counts_from_realcounts(tsim_counts = tsim_counts_mrna, 
        realcount = realcount, map_reference_cpm = model$experiment_params$map_reference_cpm, 
        map_reference_ls = model$experiment_params$map_reference_ls)
    prot_ids <- mol_info %>% filter(.data$mol == "mol_protein") %>% 
        pull(.data$val)
    count_prot_simulation <- tsim_counts[, prot_ids, drop = FALSE]
    model <- dyngen:::.add_timing(model, "7_experiment", "create output")
    sim_wcounts <- count_mrna_simulation[, model$feature_info$mol_premrna, 
        drop = FALSE]
    sim_xcounts <- count_mrna_simulation[, model$feature_info$mol_mrna, 
        drop = FALSE]
    sim_ycounts <- count_prot_simulation[, model$feature_info$mol_protein, 
        drop = FALSE]
    dimnames(sim_wcounts) <- dimnames(sim_xcounts) <- dimnames(sim_ycounts) <- list(cell_info$cell_id, 
        model$feature_info$feature_id)
    if (model$simulation_params$compute_cellwise_grn) {
        sim_cellwise_grn <- model$simulations$cellwise_grn[step_ixs, 
            , drop = FALSE]
        rownames(sim_cellwise_grn) <- cell_info$cell_id
    }
    else {
        sim_cellwise_grn <- NULL
    }
    if (model$simulation_params$compute_rna_velocity) {
        sim_rna_velocity <- model$simulations$rna_velocity[step_ixs, 
            , drop = FALSE]
        rownames(sim_rna_velocity) <- cell_info$cell_id
    }
    else {
        sim_rna_velocity <- NULL
    }
    model$experiment <- list(counts_premrna = sim_wcounts, counts_mrna = sim_xcounts, 
        counts_protein = sim_ycounts, feature_info = model$feature_info, 
        cell_info = cell_info, cellwise_grn = sim_cellwise_grn, 
        rna_velocity = sim_rna_velocity)
    dyngen:::.add_timing(model, "7_experiment", "end")
}

generate_experiment_debug_control_time = 
function (model,cell_info) 
{
    if (model$verbose) 
        cat("Simulating experiment\n")
    model <- dyngen:::.add_timing(model, "7_experiment", "sample cells")
    step_ixs <- cell_info$step_ix
    tsim_counts <- model$simulations$counts[step_ixs, , drop = FALSE]
    rownames(tsim_counts) <- cell_info$cell_id
    mol_info <- model$feature_info %>% gather("mol", "val", .data$mol_premrna, 
        .data$mol_mrna, .data$mol_protein) %>% select(.data$feature_id, 
        .data$mol, .data$val)
    model <- dyngen:::.add_timing(model, "7_experiment", "fetch realcount")
    realcount <- dyngen:::.generate_experiment_fetch_realcount(model)
    model <- dyngen:::.add_timing(model, "7_experiment", "simulate library size variation")
    mrna_ids <- mol_info %>% filter(.data$mol != "mol_protein") %>% 
        pull(.data$val)
    tsim_counts_mrna <- tsim_counts[, mrna_ids, drop = FALSE]
    count_mrna_simulation <- dyngen:::.simulate_counts_from_realcounts(tsim_counts = tsim_counts_mrna, 
        realcount = realcount, map_reference_cpm = model$experiment_params$map_reference_cpm, 
        map_reference_ls = model$experiment_params$map_reference_ls)
    prot_ids <- mol_info %>% filter(.data$mol == "mol_protein") %>% 
        pull(.data$val)
    count_prot_simulation <- tsim_counts[, prot_ids, drop = FALSE]
    model <- dyngen:::.add_timing(model, "7_experiment", "create output")
    sim_wcounts <- count_mrna_simulation[, model$feature_info$mol_premrna, 
        drop = FALSE]
    sim_xcounts <- count_mrna_simulation[, model$feature_info$mol_mrna, 
        drop = FALSE]
    sim_ycounts <- count_prot_simulation[, model$feature_info$mol_protein, 
        drop = FALSE]
    dimnames(sim_wcounts) <- dimnames(sim_xcounts) <- dimnames(sim_ycounts) <- list(cell_info$cell_id, 
        model$feature_info$feature_id)
    if (model$simulation_params$compute_cellwise_grn) {
        sim_cellwise_grn <- model$simulations$cellwise_grn[step_ixs, 
            , drop = FALSE]
        rownames(sim_cellwise_grn) <- cell_info$cell_id
    }
    else {
        sim_cellwise_grn <- NULL
    }
    if (model$simulation_params$compute_rna_velocity) {
        sim_rna_velocity <- model$simulations$rna_velocity[step_ixs, 
            , drop = FALSE]
        rownames(sim_rna_velocity) <- cell_info$cell_id
    }
    else {
        sim_rna_velocity <- NULL
    }
    model$experiment <- list(counts_premrna = sim_wcounts, counts_mrna = sim_xcounts, 
        counts_protein = sim_ycounts, feature_info = model$feature_info, 
        cell_info = cell_info, cellwise_grn = sim_cellwise_grn, 
        rna_velocity = sim_rna_velocity)
    dyngen:::.add_timing(model, "7_experiment", "end")
}

update_protein_halflife <- function(model, halflife_map) {
  feature_info <- model$feature_info
  
  for (feature_id in names(halflife_map)) {
    new_hf <- halflife_map[[feature_id]]
    
    # Update halflife
    feature_info[feature_info$feature_id == feature_id, "protein_halflife"] <- new_hf
    
    # Recompute decay rate
    decay_rate <- log(2) / new_hf
    feature_info[feature_info$feature_id == feature_id, "protein_decay_rate"] <- decay_rate
    
    # Update in simulation parameters
    param_name <- paste0("protein_decay_rate_", feature_id)
    model$simulation_system$parameters[param_name] <- as.numeric(decay_rate)
  }
  
  # Recalculate max_protein for all features
  feature_info$max_protein <- with(feature_info, translation_rate / protein_decay_rate * max_mrna)
  
  model$feature_info <- feature_info
  return(model)
}


update_adjacency_by_variance <- function(adj_matrix, counts, threshold = 0.1) {
  # Step 1: Compute variance
  variances <- apply(counts, 2, var)

  # Step 2: Filter mol_mrna_ columns with low variance
  is_mrna <- grepl("^mol_mrna_", names(variances))
  low_var_mrna <- variances[is_mrna & variances < threshold]

  # Step 3: Extract node names (remove "mol_mrna_" prefix)
  low_var_nodes <- sub("^mol_mrna_", "", names(low_var_mrna))

  # Step 4: Update adjacency matrix: set rows and columns of these nodes to 0
  common_nodes <- intersect(rownames(adj_matrix), low_var_nodes)
  adj_matrix[common_nodes, ] <- 0
  adj_matrix[, common_nodes] <- 0

  # Step 5: Return
  return(list(
    updated_adj_matrix = adj_matrix,
    filtered_nodes = common_nodes
  ))
}




evaluate_dd_graphs <- function(dd_graph, graph_truth, ground_truth_tests) {
  results_list <- list()
  
  for (time_key in names(dd_graph)) {
    time_data <- dd_graph[[time_key]]
    start_time <- as.numeric(sub("t(\\d+).*", "\\1", time_key))
    
    # Extract matrices
    inferred_raw <- time_data$inferred_raw
    inferred_adj <- time_data$inferred_adj
    raw_binary <- time_data$raw_binary
    adj_binary <- time_data$adj_binary
    
    # ----- Precision / Recall (Graph Inference) -----
    pr_raw <- compute_precision_recall_DD(inferred_raw, graph_truth)
    pr_adj <- compute_precision_recall_DD(inferred_adj, graph_truth)
    
    # ----- Test Error Rate -----
    valid_mask <- !is.na(ground_truth_tests)
    
    raw_test_error <- mean(raw_binary[valid_mask] != ground_truth_tests[valid_mask],na.rm = TRUE)
    adj_test_error <- mean(adj_binary[valid_mask] != ground_truth_tests[valid_mask],na.rm = TRUE)
    
    # ----- Collect Results -----
    results_list[[time_key]] <- data.frame(
      time_start = start_time,
      precision_raw = pr_raw$precision,
      recall_raw = pr_raw$recall,
      precision_adj = pr_adj$precision,
      recall_adj = pr_adj$recall,
      TP_raw = pr_raw$TP,
      FP_raw = pr_raw$FP,
      FN_raw = pr_raw$FN,
      TP_adj = pr_adj$TP,
      FP_adj = pr_adj$FP,
      FN_adj = pr_adj$FN,
      raw_test_error = raw_test_error,
      adj_test_error = adj_test_error
    )
  }
  
  final_df <- do.call(rbind, results_list)
  final_df <- final_df[order(final_df$time_start), ]
  rownames(final_df) <- NULL
  return(final_df)
}


compute_wilcoxon_pvals <- function(models,
                                   exp_design,
                                   wt_name,
                                   time_range,
                                   subset = FALSE,
                                   subset_threshold = 500,
                                   zero_variance_nodes = character()) {
  stopifnot(length(time_range) == 2)
  time_start <- time_range[1]
  time_end <- time_range[2]
  
  model_names <- setdiff(rownames(exp_design), wt_name)
  gene_names <- colnames(exp_design)
  
  # WT 数据选择
  model_wt <- models[[wt_name]]
  rows_selected_wt <- model_wt$simulations$meta %>%
    mutate(row_id = row_number()) %>%
    filter(sim_time >= time_start, sim_time <= time_end) %>%
    pull(row_id)
  
  if (subset && length(rows_selected_wt) > subset_threshold) {
    rows_selected_wt <- sample(rows_selected_wt, subset_threshold)
  }
  
  expr_wt <- model_wt$simulations$counts[rows_selected_wt, , drop = FALSE]
  mrna_cols <- grepl("^mol_mrna_", colnames(expr_wt))
  expr_wt <- expr_wt[, mrna_cols]
  colnames(expr_wt) <- sub("^mol_mrna_", "", colnames(expr_wt))
  
  mat_pvals <- matrix(NA, nrow = length(model_names), ncol = length(gene_names),
                      dimnames = list(model_names, gene_names))
  mat_stats <- matrix(NA, nrow = length(model_names), ncol = length(gene_names),
                      dimnames = list(model_names, gene_names))
  
  ko_indices_list <- list()
  
  for (model_name in model_names) {
    model_ko <- models[[model_name]]
    perturbed_tfs <- gene_names[exp_design[model_name, ] == 1]
    
    rows_selected_ko <- model_ko$simulations$meta %>%
      mutate(row_id = row_number()) %>%
      filter(sim_time >= time_start, sim_time <= time_end) %>%
      pull(row_id)
    
    if (subset && length(rows_selected_ko) > subset_threshold) {
      rows_selected_ko <- sample(rows_selected_ko, subset_threshold)
    }
    
    ko_indices_list[[model_name]] <- rows_selected_ko
    
    expr_ko <- model_ko$simulations$counts[rows_selected_ko, , drop = FALSE]
    expr_ko <- expr_ko[, mrna_cols]
    colnames(expr_ko) <- sub("^mol_mrna_", "", colnames(expr_ko))
    
    if (is.null(dim(expr_ko)) || is.null(dim(expr_wt))) {
      warning(paste0("Missing data for: ", model_name, " in [", time_start, ", ", time_end, "]"))
      next
    }
    
    for (gene in setdiff(gene_names, perturbed_tfs)) {
      if (!(gene %in% colnames(expr_wt)) || !(gene %in% colnames(expr_ko)) ||
          gene %in% zero_variance_nodes) {
        next
      }
      
      wt_expr <- expr_wt[, gene]
      ko_expr <- expr_ko[, gene]
      
      if (length(wt_expr) == 0 || length(ko_expr) == 0) next
      
      test_result <- tryCatch({
        wilcox.test(wt_expr, ko_expr, alternative = "two.sided")
      }, error = function(e) return(NULL))
      
      if (!is.null(test_result) && !is.na(test_result$p.value)) {
        mat_pvals[model_name, gene] <- test_result$p.value
        mat_stats[model_name, gene] <- test_result$statistic
      }
    }
  }
  
  return(list(
    pvals = mat_pvals,
    stats = mat_stats,
    wt_indices = rows_selected_wt,
    ko_indices = ko_indices_list
  ))
}

