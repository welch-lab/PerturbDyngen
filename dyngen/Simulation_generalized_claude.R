#!/usr/bin/env Rscript
# =============================================================
# Simulation_generalized_claude.R
# Claude-generated all-in-one dyngen simulation pipeline.
#
# Parameter priority (highest → lowest):
#   1. CLI --key=value arguments
#   2. --config=<path>.R  (source a cfg list, see simulation_config_claude.R)
#   3. Built-in defaults below
#
# Usage:
#   Rscript Simulation_generalized_claude.R --config=simulation_config_claude.R
#   Rscript Simulation_generalized_claude.R --config=my_cfg.R --sparsity=90
#   Rscript Simulation_generalized_claude.R --n_tfs=500 --sparsity=70 ...
#
# Parallelism:
#   Cores are read from SLURM_CPUS_PER_TASK (default 1).
#   Set #SBATCH --cpus-per-task=N in your job script to use N cores.
#
# Checkpointing:
#   Each experiment is saved individually under OUT_DIR/models_individual/.
#   On restart the script detects which files already exist and skips them,
#   so a job can be resumed without re-running completed experiments.
#
# Outputs in OUT_DIR:
#   models_all.rds            named list: "WT" + one model per experiment
#   models_individual/WT.rds  WT checkpoint
#   models_individual/<name>.rds  per-experiment checkpoints
#   EXP_design_onehot.rds    one-hot matrix  [n_exp × n_tfs]
#   EXP_design_kd.rds        KD-value matrix [n_exp × n_tfs]
#   adj_truth_binary.rds     binary adjacency  [n_tfs × n_tfs]
#   adj_truth_strength.rds   signed-weighted adjacency [n_tfs × n_tfs]
#   model_index.csv          per-experiment metadata
#   run_manifest.csv         full parameter record
# =============================================================

suppressPackageStartupMessages({
  library(dyngen)
  library(dplyr)
  library(tibble)
  library(igraph)
  library(binpackr)
  library(furrr)
  library(purrr)
  library(progress)
})

source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/dyngen_my_func.R")
source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")

# --------------- helpers -----------------------------------------------
parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]; val <- if (length(kv) >= 2) kv[2] else ""
    out[[key]] <- val
  }
  out
}

get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd()))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
}

# Coerce a string CLI value to numeric vector (comma-separated)
csv_to_num <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])

# Format a numeric KD value as a filename tag (e.g. 0.05 → "0p0500")
kd_fmt <- function(v) gsub("\\.", "p", sprintf("%.4f", v))

# Build KD-value matrix from a one-hot matrix and a per-row KD multiplier vector.
#   onehot=0 → 1  (no perturbation, full expression)
#   onehot=1 → kd_per_row[i]  (expression fraction remaining)
build_kd_matrix <- function(onehot, kd_per_row) {
  m <- matrix(1, nrow = nrow(onehot), ncol = ncol(onehot),
              dimnames = dimnames(onehot))
  for (i in seq_len(nrow(onehot))) {
    cols <- which(onehot[i, ] == 1L)
    if (length(cols) > 0L) m[i, cols] <- kd_per_row[i]
  }
  m
}

# Override cfg list with non-NULL CLI values (type-matched)
apply_cli <- function(cfg, cli) {
  for (nm in names(cli)) {
    if (!nzchar(cli[[nm]])) next
    val <- cli[[nm]]
    if (nm %in% names(cfg)) {
      orig <- cfg[[nm]]
      if (is.null(orig) || is.character(orig)) {
        cfg[[nm]] <- val
      } else if (is.logical(orig)) {
        cfg[[nm]] <- toupper(val) %in% c("TRUE", "T", "1", "YES")
      } else if (is.integer(orig)) {
        cfg[[nm]] <- as.integer(val)
      } else if (is.numeric(orig) && length(orig) == 1) {
        cfg[[nm]] <- as.numeric(val)
      } else if (is.numeric(orig) && length(orig) > 1) {
        cfg[[nm]] <- csv_to_num(val)
      }
    } else {
      cfg[[nm]] <- val   # unknown key: store as string
    }
  }
  cfg
}

# --------------- load config + apply CLI overrides --------------------
script_dir <- get_script_dir()
source(file.path(script_dir, "sampled_backbone_utils.R"))
# find_minimal_fair_ordered_K() is defined in sampled_backbone_utils.R

cli <- parse_args(commandArgs(trailingOnly = TRUE))

# Built-in defaults (mirrors simulation_config_claude.R)
cfg <- list(
  n_tfs = 500L, n_targets = 0L, n_cells = 100L, num_sims = 300L,
  sparsity = 70, force_cycle = TRUE, force_connected = FALSE,
  deadend_prob = 0.1, cycle_length_min = 3L, cycle_length_max = 5L,
  max_in = NULL, max_out = NULL, basal_value = 0.2,
  effect_probs = c(0.75, 0.25),
  strength_values = c(1, 2, 4, 10, 100),
  strength_probs  = c(0.4, 0.3, 0.15, 0.1, 0.05),
  target_states = 6L, max_reg_genes = NULL, max_terminal_states = 2L,
  reg_sign_probs = c(0.75, 0.25),
  max_nodes_exp = 8L, n_experiments = 0L,
  hub_top_k = 5L, hub_min_upstream = 3L,
  kd_effects = c(0.05, 0.1, 0.15, 0.2, 0.25),
  # kd_mode: "expand" = all KD levels × every experiment (default)
  #          "sample" = one randomly sampled KD per experiment
  kd_mode    = "expand",
  seed = 2025L, save_every = 20L, out_dir = ""
)

# Load config file if --config= provided (overrides defaults)
if (!is.null(cli$config) && nzchar(cli$config)) {
  cfg_env <- new.env(parent = baseenv())
  source(cli$config, local = cfg_env)
  if (exists("cfg", envir = cfg_env)) {
    user_cfg <- get("cfg", envir = cfg_env)
    for (nm in names(user_cfg)) cfg[[nm]] <- user_cfg[[nm]]
  }
}

# CLI overrides config (highest priority)
cfg <- apply_cli(cfg, cli)

# --------------- detect cores ------------------------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# --------------- derived values ----------------------------------------
sp       <- if (cfg$sparsity > 1) cfg$sparsity / 100 else cfg$sparsity
sp_int   <- as.integer(round(sp * 100))
sp_label <- paste0("sp", sp_int)

max_in_val  <- if (is.null(cfg$max_in))  ceiling(0.9 * cfg$n_tfs) else as.integer(cfg$max_in)
max_out_val <- if (is.null(cfg$max_out)) ceiling(0.9 * cfg$n_tfs) else as.integer(cfg$max_out)
max_reg_genes_val <- if (is.null(cfg$max_reg_genes)) ceiling(0.05 * cfg$n_tfs) else as.integer(cfg$max_reg_genes)
cycle_range <- seq(cfg$cycle_length_min, cfg$cycle_length_max)

kd_effects <- cfg$kd_effects
kd_effects <- kd_effects[!is.na(kd_effects)]
if (length(kd_effects) == 0) stop("kd_effects is empty")

OUT_DIR <- if (nzchar(cfg$out_dir)) {
  cfg$out_dir
} else {
  file.path(
    "/scratch/sigbio_project_root/sigbio_project25/wzqian/Casual_crispr/dyngen/benchmark_data",
    paste0("claude_N", cfg$n_tfs, "t", cfg$n_targets, "_", sp_label)
  )
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CACHE_DIR  <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/benchmarking/dyngen_cache/"
MODELS_DIR <- file.path(OUT_DIR, "models_individual")
dir.create(MODELS_DIR, recursive = TRUE, showWarnings = FALSE)

options(Ncpus = 8L)

message("=== Simulation_generalized_claude.R ===")
message(sprintf(
  "n_tfs=%d  n_targets=%d  n_cells=%d  num_sims=%d  %s",
  cfg$n_tfs, cfg$n_targets, cfg$n_cells, cfg$num_sims, sp_label
))
message(sprintf(
  "backbone: force_cycle=%s  force_connected=%s  deadend_prob=%.2f  cycle_range=%d:%d",
  cfg$force_cycle, cfg$force_connected, cfg$deadend_prob,
  cfg$cycle_length_min, cfg$cycle_length_max
))
message(sprintf("kd_mode=%s  kd_effects: %s", cfg$kd_mode, paste(kd_effects, collapse = ", ")))
message("out_dir:    ", OUT_DIR)
message("models_dir: ", MODELS_DIR)
message("n_cores:    ", n_cores)

# -----------------------------------------------------------------------
# STEP 1: Generate backbone + ground-truth graphs
#   Uses sample_backbone_graph() from sampled_backbone_utils.R,
#   which matches benchmark_dataset_massive.Rmd exactly:
#     target_edges = ceiling((1 - sparsity) * n_nodes * max_out)
#   with force_cycle / deadend_prob / cycle_length_range / force_connected.
# -----------------------------------------------------------------------
message("\n[1/5] Generating backbone (seed=", cfg$seed, ") ...")

backbone_list <- sample_backbone(
  n_nodes             = cfg$n_tfs,
  sparsity            = sp,
  max_in              = max_in_val,
  max_out             = max_out_val,
  force_cycle         = cfg$force_cycle,
  force_connected     = cfg$force_connected,
  deadend_prob        = cfg$deadend_prob,
  cycle_length_range  = cycle_range,
  target_states       = cfg$target_states,
  max_reg_genes       = max_reg_genes_val,
  max_terminal_states = cfg$max_terminal_states,
  basal_value         = cfg$basal_value,
  effect_probs        = cfg$effect_probs,
  strength_values     = cfg$strength_values,
  strength_probs      = cfg$strength_probs,
  reg_sign_probs      = cfg$reg_sign_probs,
  seed                = cfg$seed
)

backbone     <- backbone_list$backbone
adj_strength <- backbone_list$adj_matrix
adj_binary   <- ifelse(adj_strength != 0, 1L, 0L)

saveRDS(adj_binary,   file.path(OUT_DIR, "adj_truth_binary.rds"))
saveRDS(adj_strength, file.path(OUT_DIR, "adj_truth_strength.rds"))
message("  Observed sparsity = ", round(backbone_list$graph_stats$sparsity, 4),
        "  dead-ends = ", backbone_list$graph_stats$num_deadends,
        "  has_cycle = ", backbone_list$graph_stats$has_cycles)

# -----------------------------------------------------------------------
# STEP 2: Initialise dyngen model
# -----------------------------------------------------------------------
message("\n[2/5] Initialising dyngen model ...")

config_model <- initialise_model(
  backbone       = backbone,
  num_cells      = cfg$n_cells,
  num_targets    = cfg$n_targets,
  num_hks        = 0L,
  experiment     = experiment_synchronised(num_timepoints = 20L, pct_between = 0),
  verbose        = interactive(),
  download_cache_dir = CACHE_DIR,
  simulation_params = simulation_default(
    census_interval   = 5,
    ssa_algorithm     = ssa_etl(tau = 300 / 3600),
    experiment_params = simulation_type_wild_type(num_simulations = cfg$num_sims),
    total_time        = 400
  )
)

all_seed_common <- config_model$simulation_params$experiment_params$seed

model_common <- config_model %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard()

# scale protein half-lives x3 (standard for this pipeline)
new_half_lives <- setNames(
  model_common$feature_info$protein_halflife * 3,
  model_common$feature_info$feature_id
)
model_common <- update_protein_halflife(model_common, new_half_lives)
attr(model_common, "halflife_scaled_x3") <- TRUE

message("  Total features: ", nrow(model_common$feature_info))

# Identify TF feature IDs (used for experimental design)
tf_ids <- if ("is_tf" %in% colnames(model_common$feature_info)) {
  model_common$feature_info$feature_id[model_common$feature_info$is_tf]
} else {
  model_common$feature_info$feature_id[seq_len(cfg$n_tfs)]
}
message("  TF features for design: ", length(tf_ids))

# -----------------------------------------------------------------------
# STEP 3: Build experimental design
#   Mirrors benchmark_dataset_massive.Rmd exactly:
#     fair_ordered(N, K, 8) + build_hub_groups + multi_pass_binpackr_exp_design
#   K is auto-found (minimal) when n_experiments == 0.
# -----------------------------------------------------------------------
message("\n[3/5] Building experimental design ...")

if (cfg$n_experiments == 0L) {
  fo      <- find_minimal_fair_ordered_K(length(tf_ids), cfg$max_nodes_exp, seed = cfg$seed)
  my_exp  <- fo$result
  n_exp_K <- fo$K
} else {
  my_exp  <- fair_ordered(length(tf_ids), cfg$n_experiments, cfg$max_nodes_exp, seed = cfg$seed)
  n_exp_K <- cfg$n_experiments
}
message("  fair_ordered: K=", n_exp_K, " → ", nrow(my_exp$experiment_matrix), " random experiments")

priority_sets       <- build_hub_groups(adj_binary,
                                        max_nodes    = cfg$max_nodes_exp,
                                        top_k        = cfg$hub_top_k,
                                        min_upstream = cfg$hub_min_upstream)
priority_sets       <- unique(priority_sets)
final_priority_sets <- remove_nested_groups(priority_sets)

if (length(final_priority_sets) > 0) {
  exp_design_hub <- multi_pass_binpackr_exp_design(
    final_priority_sets,
    N_tf      = length(tf_ids),
    max_nodes = cfg$max_nodes_exp
  )
  message("  hub-group design: ", nrow(exp_design_hub), " experiments")
} else {
  exp_design_hub <- matrix(0L, nrow = 0L, ncol = length(tf_ids))
  message("  hub-group design: 0 experiments (no qualifying hubs)")
}

# Base one-hot design (n_base rows, always the same regardless of mode)
exp_design_base <- rbind(my_exp$experiment_matrix, exp_design_hub)
rownames(exp_design_base) <- c(
  paste0("EXP_",     seq_len(nrow(my_exp$experiment_matrix))),
  if (nrow(exp_design_hub) > 0L)
    paste0("HUB_EXP_", seq_len(nrow(exp_design_hub))) else character(0)
)
colnames(exp_design_base) <- tf_ids
n_base <- nrow(exp_design_base)
message("  Base design: ", n_base, " experiments x ", length(tf_ids), " TF genes")

# ---- Expand design and assign KD values by mode ----
if (cfg$kd_mode == "expand") {
  # Every base experiment is replicated once per KD level.
  # Row names encode both experiment and KD: EXP_1_KD0p0500, EXP_1_KD0p1000, ...
  n_kd     <- length(kd_effects)
  base_idx <- rep(seq_len(n_base), each = n_kd)
  kd_rep   <- rep(kd_effects, times = n_base)

  exp_design_onehot           <- exp_design_base[base_idx, , drop = FALSE]
  rownames(exp_design_onehot) <- paste0(
    rep(rownames(exp_design_base), each = n_kd), "_KD", kd_fmt(kd_rep)
  )
  kd_per_exp <- setNames(kd_rep, rownames(exp_design_onehot))
  message(sprintf("  expand mode: %d base × %d KD levels = %d total experiments",
                  n_base, n_kd, nrow(exp_design_onehot)))

} else {
  # One randomly sampled KD level per base experiment.
  # Row names: EXP_1_KD0p1000, EXP_2_KD0p0500, ...
  set.seed(cfg$seed + 1L)
  kd_samp <- sample(kd_effects, n_base, replace = TRUE)

  exp_design_onehot           <- exp_design_base
  rownames(exp_design_onehot) <- paste0(rownames(exp_design_base), "_KD", kd_fmt(kd_samp))
  kd_per_exp <- setNames(kd_samp, rownames(exp_design_onehot))
  message(sprintf("  sample mode: %d experiments, one KD sampled each", n_base))
}

n_exp_total   <- nrow(exp_design_onehot)
exp_design_kd <- build_kd_matrix(exp_design_onehot, kd_per_exp)

saveRDS(exp_design_onehot, file.path(OUT_DIR, "EXP_design_onehot.rds"))
saveRDS(exp_design_kd,     file.path(OUT_DIR, "EXP_design_kd.rds"))

# -----------------------------------------------------------------------
# Checkpoint scan: determine what still needs to be done
# -----------------------------------------------------------------------
wt_path      <- file.path(MODELS_DIR, "WT.rds")
wt_exists    <- file.exists(wt_path)

already_done <- sub("\\.rds$", "", list.files(MODELS_DIR, pattern = "\\.rds$"))
already_done <- already_done[already_done != "WT"]
todo_idx     <- which(!rownames(exp_design_onehot) %in% already_done)
n_todo       <- length(todo_idx)
n_done       <- n_exp_total - n_todo

message(sprintf(
  "\n  Checkpoint status: WT=%s  KD done=%d  KD remaining=%d",
  if (wt_exists) "cached" else "pending", n_done, n_todo
))

# -----------------------------------------------------------------------
# Precompile Gillespie reactions (needed whenever any simulation will run)
# -----------------------------------------------------------------------
need_precompile <- !wt_exists || n_todo > 0L
if (need_precompile) {
  message("\n  Precompiling Gillespie reactions (once for all experiments) ...")
  precompiled_reactions <- dyngen:::.generate_cells_precompile_reactions(model_common)
  cached_precompile     <- function(model) precompiled_reactions
  original_precompile   <- dyngen:::.generate_cells_precompile_reactions
  suppressWarnings(
    assignInNamespace(".generate_cells_precompile_reactions", cached_precompile, "dyngen")
  )
  message("  Reactions cached; DLL will be reused for all experiments.")
}

# -----------------------------------------------------------------------
# STEP 4: Wild-type simulation (with checkpoint)
# -----------------------------------------------------------------------
message("\n[4/5] Running WT simulation ...")

if (wt_exists) {
  message("  Loading WT from checkpoint ...")
  model_wt <- readRDS(wt_path)
} else {
  model_wt <- model_common
  model_wt$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = cfg$num_sims,
    timepoint       = 0.25,
    genes           = rep(list(tf_ids[seq_len(min(8L, length(tf_ids)))]), cfg$num_sims),
    num_genes       = min(8L, length(tf_ids)),
    multiplier      = 1,
    seed            = all_seed_common
  )
  model_wt <- model_wt %>% generate_cells()
  saveRDS(model_wt, wt_path)
  message("  WT done, checkpoint saved.")
}

# -----------------------------------------------------------------------
# STEP 5: KD simulations (parallel, with per-experiment checkpointing)
# Row name of exp_design_onehot is the model name (already carries KD tag).
# Each completed experiment is saved as MODELS_DIR/<model_name>.rds.
# On restart, any existing file is skipped.
# -----------------------------------------------------------------------
message("\n[5/5] Running ", n_todo, " / ", n_exp_total,
        " KD experiments (mode=", cfg$kd_mode,
        ")  cores=", n_cores, " ...")

if (n_todo > 0L) {
  # Worker: runs one KD experiment and saves its model file.
  # All variables captured from parent environment via fork.
  run_one_kd <- function(i) {
    perturbed_genes <- colnames(exp_design_onehot)[which(exp_design_onehot[i, ] == 1L)]
    if (length(perturbed_genes) == 0L) return(NULL)

    model_name <- rownames(exp_design_onehot)[i]
    kd_val     <- kd_per_exp[[model_name]]

    model_kd <- model_common
    model_kd$simulation_params$experiment_params <- simulation_type_knockdown(
      num_simulations = cfg$num_sims,
      timepoint       = 0.25,
      genes           = rep(list(perturbed_genes), cfg$num_sims),
      num_genes       = length(perturbed_genes),
      multiplier      = kd_val,
      seed            = all_seed_common
    )
    model_kd <- model_kd %>% generate_cells()

    saveRDS(model_kd, file.path(MODELS_DIR, paste0(model_name, ".rds")))
    invisible(TRUE)
  }

  results <- parallel::mclapply(
    todo_idx, run_one_kd,
    mc.cores      = n_cores,
    mc.preschedule = FALSE   # better load balance for variable-duration tasks
  )

  # Report failures
  failed_idx <- which(sapply(results, function(x) is.null(x) || inherits(x, "try-error")))
  if (length(failed_idx) > 0L) {
    warning(length(failed_idx), " experiments may have failed. Indices: ",
            paste(todo_idx[failed_idx], collapse = ", "))
  }
  message("  Parallel batch complete.")
}

# Restore original precompile function (if it was patched)
if (need_precompile) {
  suppressWarnings(
    assignInNamespace(".generate_cells_precompile_reactions", original_precompile, "dyngen")
  )
}

# -----------------------------------------------------------------------
# Assemble models_all.rds from individual checkpoint files
# -----------------------------------------------------------------------
message("\n  Assembling models_all.rds ...")

kd_names <- rownames(exp_design_onehot)
models_all <- c(
  list(WT = model_wt),
  setNames(
    lapply(kd_names, function(nm) {
      p <- file.path(MODELS_DIR, paste0(nm, ".rds"))
      if (file.exists(p)) readRDS(p) else {
        warning("Missing model file for: ", nm); NULL
      }
    }),
    kd_names
  )
)
# Drop any NULL entries (failed experiments)
models_all <- Filter(Negate(is.null), models_all)

saveRDS(models_all, file.path(OUT_DIR, "models_all.rds"))
message("  models_all.rds contains ", length(models_all), " models (WT + KD)")

# -----------------------------------------------------------------------
# Build model_index.csv from completed individual files
# -----------------------------------------------------------------------
completed_kd <- sub("\\.rds$", "", list.files(MODELS_DIR, pattern = "\\.rds$"))
completed_kd <- completed_kd[completed_kd != "WT"]

meta_df <- dplyr::bind_rows(lapply(completed_kd, function(nm) {
  if (!nm %in% rownames(exp_design_onehot)) return(NULL)
  perturbed_genes <- colnames(exp_design_onehot)[which(exp_design_onehot[nm, ] == 1L)]
  data.frame(
    model_name        = nm,
    base_experiment   = sub("_KD[0-9p]+$", "", nm),
    kd_multiplier     = kd_per_exp[[nm]],
    n_perturbed_genes = length(perturbed_genes),
    stringsAsFactors  = FALSE
  )
}))

write.csv(meta_df, file.path(OUT_DIR, "model_index.csv"), row.names = FALSE)

# Full parameter manifest
manifest <- data.frame(
  timestamp           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script              = "Simulation_generalized_claude.R",
  n_tfs               = cfg$n_tfs,
  n_targets           = cfg$n_targets,
  n_cells             = cfg$n_cells,
  num_simulations     = cfg$num_sims,
  sparsity            = sp_label,
  force_cycle         = cfg$force_cycle,
  force_connected     = cfg$force_connected,
  deadend_prob        = cfg$deadend_prob,
  cycle_length_range  = paste0(cfg$cycle_length_min, ":", cfg$cycle_length_max),
  max_in              = max_in_val,
  max_out             = max_out_val,
  basal_value         = cfg$basal_value,
  target_states       = cfg$target_states,
  max_reg_genes       = max_reg_genes_val,
  max_terminal_states = cfg$max_terminal_states,
  kd_mode             = cfg$kd_mode,
  n_fair_ordered_K    = n_exp_K,
  n_hub_experiments   = nrow(exp_design_hub),
  n_base_experiments  = n_base,
  n_total_experiments = n_exp_total,
  max_nodes_exp       = cfg$max_nodes_exp,
  hub_top_k           = cfg$hub_top_k,
  hub_min_upstream    = cfg$hub_min_upstream,
  kd_effects          = paste(kd_effects, collapse = ","),
  seed                = cfg$seed,
  n_cores             = n_cores,
  out_dir             = OUT_DIR,
  stringsAsFactors    = FALSE
)
write.csv(manifest, file.path(OUT_DIR, "run_manifest.csv"), row.names = FALSE)

message("\n=== Done ===")
message("models_all.rds         : ", file.path(OUT_DIR, "models_all.rds"))
message("models_individual/     : ", MODELS_DIR)
message("EXP_design_onehot.rds  : ", file.path(OUT_DIR, "EXP_design_onehot.rds"))
message("EXP_design_kd.rds      : ", file.path(OUT_DIR, "EXP_design_kd.rds"))
message("adj_truth_binary.rds   : ", file.path(OUT_DIR, "adj_truth_binary.rds"))
message("adj_truth_strength.rds : ", file.path(OUT_DIR, "adj_truth_strength.rds"))
message("model_index.csv        : ", file.path(OUT_DIR, "model_index.csv"))
message("run_manifest.csv       : ", file.path(OUT_DIR, "run_manifest.csv"))
