#!/usr/bin/env Rscript
# =============================================================
# Simulation_generalized_claude_perturb_seq.R
# Extension of Simulation_generalized_claude.R with a
# Perturb-seq guide RNA UMI observation layer.
#
# === Core design ===
#
# For every base experiment (perturbed gene set), the pipeline always
# runs all 5 KD levels so that guide UMI counts can be mapped to a
# meaningful KD strength via a sigmoid.
#
#   KD levels:   0.005, 0.2, 0.4, 0.6, 0.8   (multipliers into dyngen)
#   Cells/level: n_cells_per_kd (default 200)
#
# Guide UMI model  (per cell c, per target gene g, per KD level k)
# ----------------------------------------------------------------
#  1. Per-cell size factor:    s_c  ~ LogNormal(0, sigma_s^2)
#                                      sigma_s = sqrt(guide_sigma_s2)  (default sqrt(0.5))
#     E[s_c] = exp(sigma_s^2 / 2)
#
#  2. NB mean for each (cell, gene) at KD level k:
#       u_k_raw = sigmoid_inv(kd_k) / E[s_c]
#       lambda_cg = s_c * u_k_raw
#     This ensures E[X_cg] = sigmoid_inv(kd_k) so cells at each KD
#     level are centred on the target UMI for that level.
#
#     Target UMIs (from sigmoid inverse, see below):
#       kd=0.8  -> UMI_target ≈  50
#       kd=0.6  -> UMI_target ≈ 116
#       kd=0.4  -> UMI_target ≈ 170
#       kd=0.2  -> UMI_target ≈ 237
#       kd=0.005 -> UMI_target ≈ 501
#     Mean across all five levels ≈ 215.  "mean~300" refers to the
#     kd=0.2 level (237); tune guide_sigma_s2 / guide_nb_theta to taste.
#
#  3. True UMI counts:  X_cg ~ NB(lambda_cg, theta)
#
#  4. Observed UMI:     obs_cg ~ Binomial(X_cg, umi_capture_rate)
#                                  default capture = 0.90 (10% dropout)
#
# Sigmoid  (maps UMI -> continuous KD for recovery/analysis):
# ----------------------------------------------------------------
#   kd_cont(X) = 1 / (1 + exp(sigmoid_k * (X - sigmoid_u0)))
#
#   Two anchor constraints (from user spec):
#     X ≈  50 -> kd ≈ 0.8    =>  k*(50  - X0) = log(0.25) = -1.386
#     X ≈ 500 -> kd ≈ 0.005  =>  k*(500 - X0) = log(199 ) =  5.293
#   Solving: k ≈ 0.01484,  X0 ≈ 143.4   (defaults below)
#
# === Differences from base script ===
#   - kd_effects ALWAYS = c(0.005, 0.2, 0.4, 0.6, 0.8)  (fixed 5-level design)
#   - kd_mode    ALWAYS = "expand"  (all 5 levels per base experiment)
#   - n_cells refers to cells PER KD LEVEL (total/experiment = n_cells * 5)
#   - New guide UMI simulation step after dyngen runs
#
# === New CLI parameters ===
#   --n_cells_per_kd=200     cells per KD level per experiment (default 200)
#   --guide_sigma_s2=0.5     variance of log size factor (default 0.5)
#   --guide_nb_theta=10      NB dispersion for guide UMIs (default 10)
#   --sigmoid_k=0.01484      sigmoid steepness (default derived from anchors)
#   --sigmoid_u0=143.4       sigmoid midpoint  (default derived from anchors)
#   --umi_capture_rate=0.9   binomial detection probability (default 0.9)
#
# === Outputs in OUT_DIR ===
#   [inherited]
#   models_all.rds, models_individual/, EXP_design_onehot.rds,
#   EXP_design_kd.rds, adj_truth_binary.rds, adj_truth_strength.rds,
#   model_index.csv, run_manifest.csv
#
#   [new — perturb_seq/]
#   guide_params.csv           sigmoid + UMI model parameter record
#   kd_level_umi_targets.csv   per-KD-level: sigmoid-inverse UMI targets
#   guide_true_umis_all.rds    list[exp_name] -> [n_cells x n_genes] int matrix
#   guide_obs_umis_all.rds     list[exp_name] -> [n_cells x n_genes] int matrix
#   guide_kd_continuous_all.rds list[exp_name] -> [n_cells] numeric (sigmoid output)
#   cell_size_factors_all.rds  list[exp_name] -> [n_cells] numeric (s_c draws)
#   individual/                per-experiment checkpoint .rds files
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
  library(Matrix)
})

source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/dyngen_my_func.R")
source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")

# ============================================================
# SECTION 1: Helpers  (identical to base script)
# ============================================================

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) next
    kv  <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
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

csv_to_num <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])

kd_fmt <- function(v) gsub("\\.", "p", sprintf("%.4f", v))

build_kd_matrix <- function(onehot, kd_per_row) {
  m <- matrix(1, nrow = nrow(onehot), ncol = ncol(onehot),
              dimnames = dimnames(onehot))
  for (i in seq_len(nrow(onehot))) {
    cols <- which(onehot[i, ] == 1L)
    if (length(cols) > 0L) m[i, cols] <- kd_per_row[i]
  }
  m
}

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
      cfg[[nm]] <- val
    }
  }
  cfg
}

# ============================================================
# SECTION 2: Guide UMI simulation helpers  (NEW)
# ============================================================

#' Invert the sigmoid to find the UMI that maps to a given KD level.
#'   kd = 1 / (1 + exp(k * (X - X0)))  =>  X = X0 + log(1/kd - 1) / k
sigmoid_inverse <- function(kd, sigmoid_k, sigmoid_u0) {
  sigmoid_u0 + log(1 / kd - 1) / sigmoid_k
}

#' Sigmoid function: maps UMI -> continuous KD in (0, 1).
#'   High UMI -> kd -> 0 (strong knockdown)
#'   Low  UMI -> kd -> 1 (no knockdown)
sigmoid_kd <- function(umi, sigmoid_k, sigmoid_u0) {
  1 / (1 + exp(sigmoid_k * (umi - sigmoid_u0)))
}

#' Simulate guide RNA UMI counts for one (experiment x KD-level) block.
#'
#' Guide UMIs are conditioned on the true KD level of the cell so that
#' the sigmoid-recovered KD is consistent with the dyngen simulation.
#'
#' @param n_cells         Number of cells in this block.
#' @param n_target_genes  Number of perturbed genes in this experiment.
#' @param kd_level        The KD multiplier for this dyngen run (e.g. 0.005).
#' @param sigma_s         SD of log size factor  s_c ~ LogNormal(0, sigma_s^2).
#' @param nb_theta        NB dispersion (size parameter).
#' @param sigmoid_k       Sigmoid steepness.
#' @param sigmoid_u0      Sigmoid midpoint.
#' @param umi_capture_rate  Binomial detection probability per true UMI.
#' @param gene_names      Optional character vector of length n_target_genes.
#' @param cell_names      Optional character vector of length n_cells.
#' @param seed            Integer RNG seed (NULL to skip).
#'
#' @return Named list:
#'   $true_umis      [n_cells x n_genes] integer matrix  — true guide UMI counts
#'   $obs_umis       [n_cells x n_genes] integer matrix  — detected guide UMI counts
#'   $kd_continuous  [n_cells] numeric — sigmoid-inferred KD per cell
#'                   (derived from mean UMI across target genes)
#'   $s_cells        [n_cells] numeric — per-cell size factor draws
#'   $u_raw          scalar   — NB mean at unit size factor for this KD level
#'   $u_target       scalar   — sigmoid-inverse target UMI for this KD level
simulate_guide_umis <- function(
    n_cells,
    n_target_genes,
    kd_level,
    sigma_s         = sqrt(0.5),
    nb_theta        = 10.0,
    sigmoid_k       = 0.01484,
    sigmoid_u0      = 143.4,
    umi_capture_rate = 0.9,
    gene_names      = NULL,
    cell_names      = NULL,
    seed            = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (kd_level <= 0 || kd_level >= 1)
    stop("kd_level must be strictly in (0, 1) for sigmoid inversion")

  # Target mean UMI for this KD level (sigmoid inverse)
  u_target <- sigmoid_inverse(kd_level, sigmoid_k, sigmoid_u0)
  if (u_target <= 0)
    stop(sprintf("sigmoid_inverse(%.4f) = %.2f <= 0; adjust sigmoid params", kd_level, u_target))

  # Raw NB mean at unit size factor:  E[X_cg] = E[s_c] * u_raw = u_target
  E_s   <- exp(sigma_s^2 / 2)       # mean of LogNormal(0, sigma_s^2)
  u_raw <- u_target / E_s

  # Row/col labels
  row_ids  <- if (!is.null(cell_names) && length(cell_names) == n_cells) {
    cell_names
  } else paste0("cell_", seq_len(n_cells))
  col_ids  <- if (!is.null(gene_names) && length(gene_names) == n_target_genes) {
    gene_names
  } else paste0("gene_", seq_len(n_target_genes))

  # Per-cell size factors
  s_cells <- rlnorm(n_cells, meanlog = 0, sdlog = sigma_s)

  # True guide UMI counts  (per cell, per target gene, independently)
  # lambda_cg = s_c * u_raw  (same u_raw across target genes for this cell)
  true_umis <- matrix(0L, nrow = n_cells, ncol = n_target_genes,
                      dimnames = list(row_ids, col_ids))
  for (c_idx in seq_len(n_cells)) {
    lambda_c <- s_cells[c_idx] * u_raw
    true_umis[c_idx, ] <- as.integer(rnbinom(n_target_genes, size = nb_theta, mu = lambda_c))
  }

  # Observed UMI: binomial detection of each true UMI
  obs_umis <- matrix(
    as.integer(rbinom(n_cells * n_target_genes,
                      size = as.integer(true_umis),
                      prob = umi_capture_rate)),
    nrow = n_cells, ncol = n_target_genes,
    dimnames = dimnames(true_umis)
  )

  # Sigmoid-inferred KD per cell  (mean UMI across target genes -> sigmoid)
  mean_umi_per_cell <- if (n_target_genes > 1L) rowMeans(true_umis) else true_umis[, 1L]
  kd_continuous     <- sigmoid_kd(mean_umi_per_cell, sigmoid_k, sigmoid_u0)

  list(
    true_umis    = true_umis,
    obs_umis     = obs_umis,
    kd_continuous = kd_continuous,
    s_cells      = s_cells,
    u_raw        = u_raw,
    u_target     = u_target
  )
}

# ============================================================
# SECTION 3: Load config + apply CLI
# ============================================================

script_dir <- get_script_dir()
source(file.path(script_dir, "sampled_backbone_utils.R"))

cli <- parse_args(commandArgs(trailingOnly = TRUE))

cfg <- list(
  # --- backbone ---
  n_tfs = 500L, n_targets = 0L, num_sims = 300L,
  sparsity = 70, force_cycle = TRUE, force_connected = FALSE,
  deadend_prob = 0.1, cycle_length_min = 3L, cycle_length_max = 5L,
  max_in = NULL, max_out = NULL, basal_value = 0.2,
  effect_probs = c(0.75, 0.25),
  strength_values = c(1, 2, 4, 10, 100),
  strength_probs  = c(0.4, 0.3, 0.15, 0.1, 0.05),
  target_states = 6L, max_reg_genes = NULL, max_terminal_states = 2L,
  reg_sign_probs = c(0.75, 0.25),
  # --- experimental design ---
  max_nodes_exp = 8L, n_experiments = 0L,
  hub_top_k = 5L, hub_min_upstream = 3L,
  # 5 fixed KD levels — central to this script's design.
  # kd_mode is always "expand" (all 5 levels per base experiment).
  kd_effects = c(0.005, 0.2, 0.4, 0.6, 0.8),
  # n_cells_per_kd: cells per KD level per experiment (total = n_cells_per_kd * 5)
  n_cells_per_kd = 200L,
  # --- general ---
  seed = 2025L, save_every = 20L, out_dir = "",
  # --- guide UMI model (NEW) ---
  # sigma_s^2 = 0.1  =>  sigma_s = sqrt(0.1) ≈ 0.316, CV(s) ≈ 32%
  # sigma^2=0.5 (CV=81%) caused KD=0.005 cells with low s_c to fall into
  # the 150-300 UMI range (crossing the KD=0.2-0.4 territory).
  # sigma^2=0.1 keeps KD=0.005 p10 UMI ≈ 317, above the KD=0.2 mean (~237),
  # so the strongest KD level remains identifiable.
  # Note: adjacent mid-levels (0.2, 0.4, 0.6, 0.8) naturally overlap due to
  # NB dispersion regardless of sigma — this is realistic Perturb-seq behaviour.
  guide_sigma_s2   = 0.1,
  # NB dispersion; higher = less overdispersion (tighter UMI dist per level)
  guide_nb_theta   = 10.0,
  # Sigmoid anchors calibrated to user spec:
  #   UMI ≈  50 -> kd ≈ 0.8   (above 50 UMI => at least 20% KD)
  #   UMI ≈ 500 -> kd ≈ 0.005 (near 500 UMI => effect saturated)
  #   Solving: k ≈ 0.01484, X0 ≈ 143.4
  sigmoid_k        = 0.01484,
  sigmoid_u0       = 143.4,
  umi_capture_rate = 0.9
)

# Load optional config file
if (!is.null(cli$config) && nzchar(cli$config)) {
  cfg_env <- new.env(parent = baseenv())
  source(cli$config, local = cfg_env)
  if (exists("cfg", envir = cfg_env)) {
    user_cfg <- get("cfg", envir = cfg_env)
    for (nm in names(user_cfg)) cfg[[nm]] <- user_cfg[[nm]]
  }
}

cfg <- apply_cli(cfg, cli)

# Derived
sigma_s    <- sqrt(cfg$guide_sigma_s2)
E_s        <- exp(cfg$guide_sigma_s2 / 2)   # mean of LogNormal(0, sigma_s^2)

# Pre-compute per-level UMI targets (for logging and output)
kd_effects <- cfg$kd_effects
kd_effects <- kd_effects[!is.na(kd_effects)]
if (length(kd_effects) == 0) stop("kd_effects is empty")

umi_targets <- sapply(kd_effects, sigmoid_inverse, cfg$sigmoid_k, cfg$sigmoid_u0)
names(umi_targets) <- paste0("kd_", kd_effects)
if (any(umi_targets <= 0))
  stop("Some KD levels produce non-positive sigmoid-inverse UMI targets; adjust sigmoid params.")

# ============================================================
# SECTION 4: Derived values and directories
# ============================================================

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

sp       <- if (cfg$sparsity > 1) cfg$sparsity / 100 else cfg$sparsity
sp_int   <- as.integer(round(sp * 100))
sp_label <- paste0("sp", sp_int)

max_in_val        <- if (is.null(cfg$max_in))       ceiling(0.9  * cfg$n_tfs) else as.integer(cfg$max_in)
max_out_val       <- if (is.null(cfg$max_out))      ceiling(0.9  * cfg$n_tfs) else as.integer(cfg$max_out)
max_reg_genes_val <- if (is.null(cfg$max_reg_genes)) ceiling(0.05 * cfg$n_tfs) else as.integer(cfg$max_reg_genes)
cycle_range <- seq(cfg$cycle_length_min, cfg$cycle_length_max)

OUT_DIR <- if (nzchar(cfg$out_dir)) {
  cfg$out_dir
} else {
  file.path(
    "/scratch/sigbio_project_root/sigbio_project25/wzqian/Casual_crispr/dyngen/benchmark_data",
    paste0("claude_perturb_seq_N", cfg$n_tfs, "t", cfg$n_targets, "_", sp_label)
  )
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CACHE_DIR       <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/benchmarking/dyngen_cache/"
MODELS_DIR      <- file.path(OUT_DIR, "models_individual")
PERTURB_DIR     <- file.path(OUT_DIR, "perturb_seq")
PERTURB_IND_DIR <- file.path(PERTURB_DIR, "individual")

dir.create(MODELS_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(PERTURB_IND_DIR, recursive = TRUE, showWarnings = FALSE)

options(Ncpus = 8L)

message("=== Simulation_generalized_claude_perturb_seq.R ===")
message(sprintf(
  "n_tfs=%d  n_targets=%d  n_cells_per_kd=%d  num_sims=%d  %s",
  cfg$n_tfs, cfg$n_targets, cfg$n_cells_per_kd, cfg$num_sims, sp_label
))
message(sprintf(
  "KD levels (%d): %s  [always expand, all 5 per experiment]",
  length(kd_effects), paste(kd_effects, collapse = ", ")
))
message("  Sigmoid:  kd(X) = 1 / (1 + exp(k*(X-X0)))")
message(sprintf("    k=%.5f  X0=%.2f", cfg$sigmoid_k, cfg$sigmoid_u0))
message("  UMI targets per KD level (sigmoid inverse):")
for (i in seq_along(kd_effects)) {
  message(sprintf("    kd=%-6s -> UMI_target=%.1f  u_raw=%.1f",
                  kd_effects[i], umi_targets[i], umi_targets[i] / E_s))
}
message(sprintf(
  "  Guide model: sigma_s=%.3f (sigma_s2=%.2f)  nb_theta=%.1f  capture_rate=%.2f",
  sigma_s, cfg$guide_sigma_s2, cfg$guide_nb_theta, cfg$umi_capture_rate
))
message("out_dir:    ", OUT_DIR)
message("n_cores:    ", n_cores)

# ============================================================
# STEP 1: Backbone + ground-truth graphs
# ============================================================
message("\n[1/6] Generating backbone (seed=", cfg$seed, ") ...")

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
message("  Sparsity=", round(backbone_list$graph_stats$sparsity, 4),
        "  dead-ends=", backbone_list$graph_stats$num_deadends,
        "  has_cycle=", backbone_list$graph_stats$has_cycles)

# ============================================================
# STEP 2: Initialise dyngen model
# ============================================================
message("\n[2/6] Initialising dyngen model ...")

# n_cells for dyngen = n_cells_per_kd (one run per KD level)
config_model <- initialise_model(
  backbone       = backbone,
  num_cells      = cfg$n_cells_per_kd,
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

new_half_lives <- setNames(
  model_common$feature_info$protein_halflife * 3,
  model_common$feature_info$feature_id
)
model_common <- update_protein_halflife(model_common, new_half_lives)
attr(model_common, "halflife_scaled_x3") <- TRUE

tf_ids <- if ("is_tf" %in% colnames(model_common$feature_info)) {
  model_common$feature_info$feature_id[model_common$feature_info$is_tf]
} else {
  model_common$feature_info$feature_id[seq_len(cfg$n_tfs)]
}
message("  Total features=", nrow(model_common$feature_info),
        "  TF features=", length(tf_ids))

# ============================================================
# STEP 3: Experimental design
#   Always: expand all 5 KD levels per base experiment.
#   Model name: {BASE_EXP_NAME}_KD{kd_fmt}
# ============================================================
message("\n[3/6] Building experimental design ...")

if (cfg$n_experiments == 0L) {
  fo      <- find_minimal_fair_ordered_K(length(tf_ids), cfg$max_nodes_exp, seed = cfg$seed)
  my_exp  <- fo$result
  n_exp_K <- fo$K
} else {
  my_exp  <- fair_ordered(length(tf_ids), cfg$n_experiments, cfg$max_nodes_exp, seed = cfg$seed)
  n_exp_K <- cfg$n_experiments
}
message("  fair_ordered: K=", n_exp_K,
        "  random experiments=", nrow(my_exp$experiment_matrix))

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

exp_design_base <- rbind(my_exp$experiment_matrix, exp_design_hub)
rownames(exp_design_base) <- c(
  paste0("EXP_",     seq_len(nrow(my_exp$experiment_matrix))),
  if (nrow(exp_design_hub) > 0L)
    paste0("HUB_EXP_", seq_len(nrow(exp_design_hub))) else character(0)
)
colnames(exp_design_base) <- tf_ids
n_base <- nrow(exp_design_base)

# Expand: every base experiment x every KD level
n_kd     <- length(kd_effects)
base_idx <- rep(seq_len(n_base), each = n_kd)
kd_rep   <- rep(kd_effects, times = n_base)

exp_design_onehot           <- exp_design_base[base_idx, , drop = FALSE]
rownames(exp_design_onehot) <- paste0(
  rep(rownames(exp_design_base), each = n_kd), "_KD", kd_fmt(kd_rep)
)
kd_per_exp    <- setNames(kd_rep,            rownames(exp_design_onehot))
n_exp_total   <- nrow(exp_design_onehot)
exp_design_kd <- build_kd_matrix(exp_design_onehot, kd_per_exp)

message(sprintf("  Total design: %d base x %d KD levels = %d experiments",
                n_base, n_kd, n_exp_total))

saveRDS(exp_design_onehot, file.path(OUT_DIR, "EXP_design_onehot.rds"))
saveRDS(exp_design_kd,     file.path(OUT_DIR, "EXP_design_kd.rds"))

# ============================================================
# Checkpoint scan
# ============================================================
wt_path   <- file.path(MODELS_DIR, "WT.rds")
wt_exists <- file.exists(wt_path)

already_done <- sub("\\.rds$", "", list.files(MODELS_DIR, pattern = "\\.rds$"))
already_done <- already_done[already_done != "WT"]
todo_idx     <- which(!rownames(exp_design_onehot) %in% already_done)
n_todo       <- length(todo_idx)

already_perturbed <- sub("\\.rds$", "", list.files(PERTURB_IND_DIR, pattern = "\\.rds$"))
todo_perturb_idx  <- which(!rownames(exp_design_onehot) %in% already_perturbed)

message(sprintf(
  "\n  Checkpoint: WT=%s  dyngen done=%d / %d  guide done=%d / %d",
  if (wt_exists) "cached" else "pending",
  n_exp_total - n_todo, n_exp_total,
  n_exp_total - length(todo_perturb_idx), n_exp_total
))

# ============================================================
# Precompile Gillespie reactions
# ============================================================
need_precompile <- !wt_exists || n_todo > 0L
if (need_precompile) {
  message("\n  Precompiling Gillespie reactions ...")
  precompiled_reactions <- dyngen:::.generate_cells_precompile_reactions(model_common)
  cached_precompile     <- function(model) precompiled_reactions
  original_precompile   <- dyngen:::.generate_cells_precompile_reactions
  suppressWarnings(
    assignInNamespace(".generate_cells_precompile_reactions", cached_precompile, "dyngen")
  )
}

# ============================================================
# STEP 4: Wild-type simulation
# ============================================================
message("\n[4/6] Running WT simulation ...")

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
  message("  WT done.")
}

# ============================================================
# STEP 5: KD dyngen simulations
#   KD multiplier passed to dyngen = kd_per_exp[[model_name]]
#   This IS the final effective KD (no guide-factor scaling needed).
# ============================================================
message("\n[5/6] Running ", n_todo, " / ", n_exp_total,
        " KD simulations  cores=", n_cores, " ...")

if (n_todo > 0L) {
  run_one_kd <- function(i) {
    model_name      <- rownames(exp_design_onehot)[i]
    perturbed_genes <- colnames(exp_design_onehot)[which(exp_design_onehot[i, ] == 1L)]
    if (length(perturbed_genes) == 0L) return(NULL)

    kd_val   <- kd_per_exp[[model_name]]
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
    mc.cores       = n_cores,
    mc.preschedule = FALSE
  )

  failed_idx <- which(sapply(results, function(x) is.null(x) || inherits(x, "try-error")))
  if (length(failed_idx) > 0L)
    warning(length(failed_idx), " experiments failed. Indices: ",
            paste(todo_idx[failed_idx], collapse = ", "))
  message("  KD batch complete.")
}

if (need_precompile) {
  suppressWarnings(
    assignInNamespace(".generate_cells_precompile_reactions", original_precompile, "dyngen")
  )
}

# ============================================================
# STEP 6: Guide RNA UMI simulation  (NEW)
#
# For each completed (base_exp x KD_level) block:
#   - Identify perturbed target genes (n_genes)
#   - Identify number of cells in the dyngen model (n_cells_per_kd)
#   - Call simulate_guide_umis() with kd_level-specific u_raw
#   - Save individual checkpoint to PERTURB_IND_DIR/
#
# The per-cell sigmoid-inferred kd_continuous is saved alongside
# true and observed UMI matrices so downstream analyses can test
# how well guide UMI counts recover the true KD label.
# ============================================================
message("\n[6/6] Simulating guide RNA UMI counts (",
        length(todo_perturb_idx), " / ", n_exp_total, " remaining) ...")

get_cell_names_from_model <- function(model) {
  if (!is.null(model$counts) && !is.null(rownames(model$counts)))
    return(rownames(model$counts))
  if (!is.null(model$cell_info) && !is.null(model$cell_info$cell_id))
    return(model$cell_info$cell_id)
  paste0("cell_", seq_len(cfg$n_cells_per_kd))
}

run_one_guide <- function(i) {
  model_name      <- rownames(exp_design_onehot)[i]
  perturbed_genes <- colnames(exp_design_onehot)[which(exp_design_onehot[i, ] == 1L)]
  kd_val          <- kd_per_exp[[model_name]]

  if (length(perturbed_genes) == 0L) return(NULL)

  model_path <- file.path(MODELS_DIR, paste0(model_name, ".rds"))
  if (!file.exists(model_path)) {
    warning("Missing dyngen model for: ", model_name, " — skipping guide UMI sim")
    return(NULL)
  }

  model_kd   <- readRDS(model_path)
  cell_names <- get_cell_names_from_model(model_kd)
  n_cells_i  <- length(cell_names)

  # Unique seed per experiment
  exp_seed <- cfg$seed + i * 1000L

  ps <- simulate_guide_umis(
    n_cells          = n_cells_i,
    n_target_genes   = length(perturbed_genes),
    kd_level         = kd_val,
    sigma_s          = sigma_s,
    nb_theta         = cfg$guide_nb_theta,
    sigmoid_k        = cfg$sigmoid_k,
    sigmoid_u0       = cfg$sigmoid_u0,
    umi_capture_rate = cfg$umi_capture_rate,
    gene_names       = perturbed_genes,
    cell_names       = cell_names,
    seed             = exp_seed
  )

  result <- list(
    model_name      = model_name,
    perturbed_genes = perturbed_genes,
    kd_level        = kd_val,
    n_cells         = n_cells_i,
    u_target        = ps$u_target,    # sigmoid-inverse UMI for this KD level
    u_raw           = ps$u_raw,       # NB mean at unit size factor
    true_umis       = ps$true_umis,   # [n_cells x n_genes]
    obs_umis        = ps$obs_umis,    # [n_cells x n_genes]
    kd_continuous   = ps$kd_continuous, # [n_cells] sigmoid-inferred KD
    s_cells         = ps$s_cells      # [n_cells] size factor draws
  )

  saveRDS(result, file.path(PERTURB_IND_DIR, paste0(model_name, ".rds")))
  invisible(TRUE)
}

for (i in todo_perturb_idx) {
  nm <- rownames(exp_design_onehot)[i]
  tryCatch(
    run_one_guide(i),
    error = function(e) warning("Guide UMI sim failed for ", nm, ": ", e$message)
  )
}
message("  Guide UMI simulation complete.")

# ============================================================
# Assemble models_all.rds
# ============================================================
message("\n  Assembling models_all.rds ...")

kd_names <- rownames(exp_design_onehot)
models_all <- c(
  list(WT = model_wt),
  setNames(
    lapply(kd_names, function(nm) {
      p <- file.path(MODELS_DIR, paste0(nm, ".rds"))
      if (file.exists(p)) readRDS(p) else { warning("Missing: ", nm); NULL }
    }),
    kd_names
  )
)
models_all <- Filter(Negate(is.null), models_all)
saveRDS(models_all, file.path(OUT_DIR, "models_all.rds"))
message("  models_all.rds: ", length(models_all), " entries (WT + KD)")

# ============================================================
# Assemble guide UMI combined output files
# ============================================================
message("\n  Assembling perturb_seq combined outputs ...")

completed_kd <- sub("\\.rds$", "", list.files(MODELS_DIR, pattern = "\\.rds$"))
completed_kd <- completed_kd[completed_kd != "WT"]

perturb_results <- setNames(
  lapply(completed_kd, function(nm) {
    p <- file.path(PERTURB_IND_DIR, paste0(nm, ".rds"))
    if (file.exists(p)) readRDS(p) else NULL
  }),
  completed_kd
)
perturb_results <- Filter(Negate(is.null), perturb_results)
n_perturb_done  <- length(perturb_results)
message("  Guide UMI results: ", n_perturb_done, " / ", length(completed_kd), " experiments")

saveRDS(lapply(perturb_results, `[[`, "true_umis"),
        file.path(PERTURB_DIR, "guide_true_umis_all.rds"))
saveRDS(lapply(perturb_results, `[[`, "obs_umis"),
        file.path(PERTURB_DIR, "guide_obs_umis_all.rds"))
saveRDS(lapply(perturb_results, `[[`, "kd_continuous"),
        file.path(PERTURB_DIR, "guide_kd_continuous_all.rds"))
saveRDS(lapply(perturb_results, `[[`, "s_cells"),
        file.path(PERTURB_DIR, "cell_size_factors_all.rds"))

# ============================================================
# KD level UMI target table
# ============================================================
kd_umi_table <- data.frame(
  kd_level   = kd_effects,
  u_target   = umi_targets,
  u_raw      = umi_targets / E_s,
  description = c(
    "Near-complete KD: UMI~501 (saturated region of sigmoid)",
    "Strong KD: UMI~237",
    "Moderate KD: UMI~170",
    "Partial KD: UMI~116",
    "Weak KD: UMI~50 (above 50 => at least 20% KD)"
  )[order(kd_effects)],
  stringsAsFactors = FALSE
)
write.csv(kd_umi_table, file.path(PERTURB_DIR, "kd_level_umi_targets.csv"), row.names = FALSE)

# ============================================================
# Guide model parameter record
# ============================================================
guide_params <- data.frame(
  parameter   = c("kd_levels", "n_cells_per_kd",
                  "guide_sigma_s2", "guide_nb_theta",
                  "sigmoid_k", "sigmoid_u0",
                  "umi_capture_rate",
                  "size_factor_model", "umi_count_model", "detection_model"),
  value       = c(paste(kd_effects, collapse = ","),
                  cfg$n_cells_per_kd,
                  cfg$guide_sigma_s2, cfg$guide_nb_theta,
                  cfg$sigmoid_k, cfg$sigmoid_u0,
                  cfg$umi_capture_rate,
                  sprintf("LogNormal(0, %.2f)", cfg$guide_sigma_s2),
                  sprintf("NB(mu = s_c * u_raw, size = %.1f)", cfg$guide_nb_theta),
                  sprintf("Binomial(true_umi, p=%.2f)", cfg$umi_capture_rate)),
  description = c(
    "5 fixed KD multipliers passed to dyngen",
    "Cells per KD level per experiment (total = n_cells_per_kd * 5)",
    "Variance of log size factor: s_c ~ LogNormal(0, sigma_s2)",
    "NB dispersion (higher = tighter UMI distribution per level)",
    "Sigmoid steepness: kd(X) = 1/(1+exp(k*(X-X0)))",
    "Sigmoid midpoint (UMI value giving kd=0.5)",
    "Binomial detection probability per true UMI (90% capture)",
    "Per-cell guide uptake scaling",
    "True guide UMI counts per cell per target gene",
    "Observed guide UMI counts after droplet capture"
  ),
  stringsAsFactors = FALSE
)
write.csv(guide_params, file.path(PERTURB_DIR, "guide_params.csv"), row.names = FALSE)

# ============================================================
# model_index.csv
# ============================================================
completed_valid <- completed_kd[completed_kd %in% rownames(exp_design_onehot)]
meta_df <- dplyr::bind_rows(lapply(completed_valid, function(nm) {
  perturbed_genes <- colnames(exp_design_onehot)[which(exp_design_onehot[nm, ] == 1L)]
  has_guide       <- nm %in% names(perturb_results)
  kd_val          <- kd_per_exp[[nm]]
  data.frame(
    model_name        = nm,
    base_experiment   = sub("_KD[0-9p]+$", "", nm),
    kd_level          = kd_val,
    umi_target        = sigmoid_inverse(kd_val, cfg$sigmoid_k, cfg$sigmoid_u0),
    n_perturbed_genes = length(perturbed_genes),
    n_cells           = cfg$n_cells_per_kd,
    has_guide_umis    = has_guide,
    stringsAsFactors  = FALSE
  )
}))
write.csv(meta_df, file.path(OUT_DIR, "model_index.csv"), row.names = FALSE)

# ============================================================
# run_manifest.csv
# ============================================================
manifest <- data.frame(
  timestamp            = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script               = "Simulation_generalized_claude_perturb_seq.R",
  n_tfs                = cfg$n_tfs,
  n_targets            = cfg$n_targets,
  n_cells_per_kd       = cfg$n_cells_per_kd,
  n_cells_total_per_exp = cfg$n_cells_per_kd * length(kd_effects),
  num_simulations      = cfg$num_sims,
  sparsity             = sp_label,
  force_cycle          = cfg$force_cycle,
  force_connected      = cfg$force_connected,
  deadend_prob         = cfg$deadend_prob,
  cycle_length_range   = paste0(cfg$cycle_length_min, ":", cfg$cycle_length_max),
  max_in               = max_in_val,
  max_out              = max_out_val,
  basal_value          = cfg$basal_value,
  target_states        = cfg$target_states,
  max_reg_genes        = max_reg_genes_val,
  max_terminal_states  = cfg$max_terminal_states,
  kd_effects           = paste(kd_effects, collapse = ","),
  n_base_experiments   = n_base,
  n_total_experiments  = n_exp_total,
  max_nodes_exp        = cfg$max_nodes_exp,
  hub_top_k            = cfg$hub_top_k,
  hub_min_upstream     = cfg$hub_min_upstream,
  guide_sigma_s2       = cfg$guide_sigma_s2,
  guide_nb_theta       = cfg$guide_nb_theta,
  sigmoid_k            = cfg$sigmoid_k,
  sigmoid_u0           = cfg$sigmoid_u0,
  umi_capture_rate     = cfg$umi_capture_rate,
  seed                 = cfg$seed,
  n_cores              = n_cores,
  out_dir              = OUT_DIR,
  n_guide_complete     = n_perturb_done,
  stringsAsFactors     = FALSE
)
write.csv(manifest, file.path(OUT_DIR, "run_manifest.csv"), row.names = FALSE)

# ============================================================
# Final summary
# ============================================================
message("\n=== Done ===")
message("models_all.rds             : ", file.path(OUT_DIR, "models_all.rds"),
        "  (", length(models_all), " models)")
message("EXP_design_onehot.rds      : ", file.path(OUT_DIR, "EXP_design_onehot.rds"))
message("EXP_design_kd.rds          : ", file.path(OUT_DIR, "EXP_design_kd.rds"))
message("adj_truth_{binary,strength}: ", OUT_DIR)
message("model_index.csv            : ", file.path(OUT_DIR, "model_index.csv"))
message("run_manifest.csv           : ", file.path(OUT_DIR, "run_manifest.csv"))
message("--- perturb_seq/ ---")
message("guide_params.csv           : ", file.path(PERTURB_DIR, "guide_params.csv"))
message("kd_level_umi_targets.csv   : ", file.path(PERTURB_DIR, "kd_level_umi_targets.csv"))
message("guide_true_umis_all.rds    : ", file.path(PERTURB_DIR, "guide_true_umis_all.rds"))
message("guide_obs_umis_all.rds     : ", file.path(PERTURB_DIR, "guide_obs_umis_all.rds"))
message("guide_kd_continuous_all    : ", file.path(PERTURB_DIR, "guide_kd_continuous_all.rds"))
message("cell_size_factors_all      : ", file.path(PERTURB_DIR, "cell_size_factors_all.rds"))
message("individual/                : ", PERTURB_IND_DIR,
        "  (", n_perturb_done, " experiments)")
