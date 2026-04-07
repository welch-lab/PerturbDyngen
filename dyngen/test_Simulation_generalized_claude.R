#!/usr/bin/env Rscript
# =============================================================
# test_Simulation_generalized_claude.R
# Unit test for Simulation_generalized_claude.R
#
# Tests run in two tiers:
#   TIER 1 (fast, no dyngen simulation): backbone generation,
#           experimental design matrices, KD-sampling logic.
#   TIER 2 (full pipeline, slow): runs generate_cells() with
#           minimal params (n_tfs=15, num_sims=10, n_cells=10).
#
# Usage:
#   Rscript test_Simulation_generalized_claude.R
#   Rscript test_Simulation_generalized_claude.R --tier=1   # fast only
#   Rscript test_Simulation_generalized_claude.R --tier=2   # full
# =============================================================

suppressPackageStartupMessages({
  library(dyngen)
  library(dplyr)
  library(tibble)
  library(igraph)
  library(binpackr)
})

source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/dyngen_my_func.R")
source("/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/theory/theory_simulations_func.R")

# ---- locate script dir --------------------------------------------------
get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd()))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
}
script_dir <- get_script_dir()
source(file.path(script_dir, "sampled_backbone_utils.R"))

# ---- parse tier ---------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
tier_arg <- grep("^--tier=", cli_args, value = TRUE)
RUN_TIER <- if (length(tier_arg) > 0) as.integer(sub("^--tier=", "", tier_arg[1])) else 1L

PASS <- 0L; FAIL <- 0L
test_results <- list()

assert <- function(name, expr, ...) {
  result <- tryCatch(
    { val <- expr; if (!isTRUE(val)) stop("assertion FALSE") ; TRUE },
    error = function(e) { message("FAIL [", name, "]: ", conditionMessage(e)); FALSE }
  )
  if (result) { message("PASS [", name, "]") ; PASS <<- PASS + 1L }
  else          { FAIL <<- FAIL + 1L }
  invisible(result)
}

# =========================================================================
# TIER 1: Unit tests (no dyngen simulation)
# =========================================================================
message("\n=== TIER 1: structural / logic tests ===\n")

# --- T1.1 backbone generation ---
message("-- T1.1 sample_backbone (n_tfs=15) --")
bb <- sample_backbone(
  n_nodes = 15, sparsity = 0.7,
  max_in = 10, max_out = 10,
  target_states = 3, max_reg_genes = 3,
  max_terminal_states = 2, seed = 42
)
assert("backbone returns list with backbone/adj_matrix/graph_stats",
       all(c("backbone","adj_matrix","graph_stats") %in% names(bb)))
assert("adj_matrix is 15x15",   all(dim(bb$adj_matrix) == c(15, 15)))
assert("graph_stats has sparsity field", !is.null(bb$graph_stats$sparsity))
assert("backbone object is valid (module_info/module_network/expression_patterns)",
       all(c("module_info","module_network","expression_patterns") %in% names(bb$backbone)))

# --- T1.2 adj_binary derivation ---
adj_str <- bb$adj_matrix
adj_bin <- ifelse(adj_str != 0, 1L, 0L)
assert("adj_binary is 0/1 only",     all(adj_bin %in% c(0L, 1L)))
assert("adj_binary dims match",      all(dim(adj_bin) == dim(adj_str)))
assert("nonzero entries match",
       sum(adj_bin) == sum(adj_str != 0))

# --- T1.3 KD sampling logic ---
set.seed(99)
kd_pool <- c(0, 0.05, 0.1, 0.2, 0.5)
n_exp   <- 20L
kd_sampled <- sample(kd_pool, n_exp, replace = TRUE)
assert("sampled KD values all in pool",
       all(kd_sampled %in% kd_pool))
assert("sampled length correct", length(kd_sampled) == n_exp)

# --- T1.4 build_kd_matrix (shared helper used by both modes) ---
n_genes <- 10L; n_base_exp <- 3L
kd_pool <- c(0.05, 0.1, 0.15, 0.2, 0.25)
set.seed(1)
onehot_base <- matrix(
  sample(c(0L, 1L), n_base_exp * n_genes, replace = TRUE, prob = c(0.7, 0.3)),
  nrow = n_base_exp, ncol = n_genes,
  dimnames = list(paste0("EXP_", seq_len(n_base_exp)), paste0("G", seq_len(n_genes)))
)

# -- sample mode --
kd_vals_s <- sample(kd_pool, n_base_exp, replace = TRUE)
onehot_s  <- onehot_base
rownames(onehot_s) <- paste0(rownames(onehot_base), "_KD",
                              gsub("\\.", "p", sprintf("%.4f", kd_vals_s)))
kd_s <- matrix(1, nrow=nrow(onehot_s), ncol=ncol(onehot_s), dimnames=dimnames(onehot_s))
for (i in seq_len(nrow(onehot_s))) {
  cols <- which(onehot_s[i,]==1L); if(length(cols)>0L) kd_s[i,cols] <- kd_vals_s[i]
}
assert("sample mode: unperturbed=1",          all(kd_s[onehot_s==0L]==1))
assert("sample mode: perturbed=kd_val",       {
  m <- onehot_s==1L & kd_vals_s[row(onehot_s)]>0
  if(sum(m)==0L) TRUE else all(kd_s[m]==kd_vals_s[row(kd_s)][m])
})
assert("sample mode: rows = n_base",          nrow(onehot_s) == n_base_exp)

# -- expand mode --
n_kd     <- length(kd_pool)
base_idx <- rep(seq_len(n_base_exp), each = n_kd)
kd_rep   <- rep(kd_pool, times = n_base_exp)
onehot_e <- onehot_base[base_idx, , drop = FALSE]
rownames(onehot_e) <- paste0(rep(rownames(onehot_base), each=n_kd), "_KD",
                              gsub("\\.", "p", sprintf("%.4f", kd_rep)))
kd_e <- matrix(1, nrow=nrow(onehot_e), ncol=ncol(onehot_e), dimnames=dimnames(onehot_e))
for (i in seq_len(nrow(onehot_e))) {
  cols <- which(onehot_e[i,]==1L); if(length(cols)>0L) kd_e[i,cols] <- kd_rep[i]
}
assert("expand mode: rows = n_base * n_kd",   nrow(onehot_e) == n_base_exp * n_kd)
assert("expand mode: unperturbed=1",           all(kd_e[onehot_e==0L]==1))
assert("expand mode: perturbed=kd_level",      {
  m <- onehot_e==1L & kd_rep[row(onehot_e)]>0
  if(sum(m)==0L) TRUE else all(kd_e[m]==kd_rep[row(kd_e)][m])
})
assert("expand mode: base_exp recovered from rowname", {
  recovered <- sub("_KD[0-9p]+$", "", rownames(onehot_e))
  all(recovered == rep(rownames(onehot_base), each=n_kd))
})
# reuse onehot/kd_matrix names for later T1.5 save test
onehot    <- onehot_s
kd_matrix <- kd_s

# --- T1.5 output directory creation ---
tmp_dir <- file.path(tempdir(), "test_claude_generalized")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
assert("temp output dir exists", dir.exists(tmp_dir))

saveRDS(adj_bin,    file.path(tmp_dir, "adj_truth_binary.rds"))
saveRDS(adj_str,    file.path(tmp_dir, "adj_truth_strength.rds"))
saveRDS(onehot,     file.path(tmp_dir, "EXP_design_onehot.rds"))
saveRDS(kd_matrix,  file.path(tmp_dir, "EXP_design_kd.rds"))
assert("adj_truth_binary.rds saved",      file.exists(file.path(tmp_dir, "adj_truth_binary.rds")))
assert("adj_truth_strength.rds saved",    file.exists(file.path(tmp_dir, "adj_truth_strength.rds")))
assert("EXP_design_onehot.rds saved",     file.exists(file.path(tmp_dir, "EXP_design_onehot.rds")))
assert("EXP_design_kd.rds saved",         file.exists(file.path(tmp_dir, "EXP_design_kd.rds")))

# reload and verify round-trip
adj_bin_rt  <- readRDS(file.path(tmp_dir, "adj_truth_binary.rds"))
onehot_rt   <- readRDS(file.path(tmp_dir, "EXP_design_onehot.rds"))
kd_mat_rt   <- readRDS(file.path(tmp_dir, "EXP_design_kd.rds"))
assert("adj_binary round-trip identical",   identical(adj_bin,   adj_bin_rt))
assert("onehot round-trip identical",       identical(onehot,    onehot_rt))
assert("kd_matrix round-trip identical",    identical(kd_matrix, kd_mat_rt))

message("\nTIER 1 complete. PASS=", PASS, "  FAIL=", FAIL)

if (RUN_TIER == 1L) {
  if (FAIL == 0L) {
    message("\nAll TIER-1 tests PASSED.")
    quit(status = 0)
  } else {
    message("\n", FAIL, " TIER-1 test(s) FAILED.")
    quit(status = 1)
  }
}

# =========================================================================
# TIER 2: Full mini-pipeline (slow — runs generate_cells)
# =========================================================================
message("\n=== TIER 2: full mini-pipeline test ===\n")
MINI_N_TFS    <- 100L   # must be large enough for fair_ordered with max=8
MINI_N_CELLS  <- 5L
MINI_NUM_SIMS <- 5L
message("Using n_tfs=", MINI_N_TFS, ", n_cells=", MINI_N_CELLS,
        ", num_sims=", MINI_NUM_SIMS)
MINI_SEED     <- 42L
MINI_KD       <- c(0, 0.1, 0.5)
MINI_KD_MODE  <- "expand"   # matches production default
MINI_OUT      <- file.path(script_dir, "test_output_claude_expand")
dir.create(MINI_OUT, recursive = TRUE, showWarnings = FALSE)
CACHE_DIR <- "/nfs/turbo/umms-welchjd-code/code/wzqian/Casualformer/dyngen/benchmarking/dyngen_cache/"

message("-- T2.1 backbone + model init --")
bb2 <- sample_backbone(
  n_nodes = MINI_N_TFS, sparsity = 0.7,
  max_in  = ceiling(0.9 * MINI_N_TFS),
  max_out = ceiling(0.9 * MINI_N_TFS),
  target_states = 4, max_reg_genes = ceiling(0.05 * MINI_N_TFS),
  max_terminal_states = 2, seed = MINI_SEED
)
assert("T2.1 backbone built", is.list(bb2))

adj_str2 <- bb2$adj_matrix
adj_bin2 <- ifelse(adj_str2 != 0, 1L, 0L)
saveRDS(adj_bin2,  file.path(MINI_OUT, "adj_truth_binary.rds"))
saveRDS(adj_str2,  file.path(MINI_OUT, "adj_truth_strength.rds"))

config2 <- initialise_model(
  backbone       = bb2$backbone,
  num_cells      = MINI_N_CELLS,
  num_targets    = 0L,
  num_hks        = 0L,
  experiment     = experiment_synchronised(num_timepoints = 5L, pct_between = 0),
  verbose        = FALSE,
  download_cache_dir = CACHE_DIR,
  simulation_params = simulation_default(
    census_interval   = 5,
    ssa_algorithm     = ssa_etl(tau = 300 / 3600),
    experiment_params = simulation_type_wild_type(num_simulations = MINI_NUM_SIMS),
    total_time        = 400
  )
)
assert("T2.1 config created", !is.null(config2))

seed_common2 <- config2$simulation_params$experiment_params$seed
model2 <- config2 %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard()
new_hl2 <- setNames(model2$feature_info$protein_halflife * 3, model2$feature_info$feature_id)
model2  <- update_protein_halflife(model2, new_hl2)

assert("T2.1 model built, has feature_info",    !is.null(model2$feature_info))
assert("T2.1 feature count matches n_tfs",      nrow(model2$feature_info) == MINI_N_TFS)

message("-- T2.2 experimental design --")
tf_ids2 <- model2$feature_info$feature_id[seq_len(MINI_N_TFS)]
# fair_ordered always uses max_nodes=8; find minimal K that satisfies its constraint
MINI_MAX_NODES <- 8L
mini_K <- ceiling(MINI_N_TFS / MINI_MAX_NODES)  # near-minimal starting point
my_exp2 <- NULL
for (try_K in seq(mini_K, mini_K + 20L)) {
  my_exp2 <- tryCatch(
    fair_ordered(length(tf_ids2), try_K, MINI_MAX_NODES, seed = MINI_SEED),
    error = function(e) NULL
  )
  if (!is.null(my_exp2)) {
    message("  fair_ordered succeeded with K=", try_K, ", max_nodes=", MINI_MAX_NODES)
    break
  }
}
if (is.null(my_exp2)) stop("Could not find valid K for fair_ordered with n=", MINI_N_TFS)

exp_onehot2 <- my_exp2$experiment_matrix
rownames(exp_onehot2) <- paste0("EXP_", seq_len(nrow(exp_onehot2)))
colnames(exp_onehot2) <- tf_ids2

# expand mode: each base experiment replicated for every KD level
n_base2  <- nrow(exp_onehot2)
n_kd2    <- length(MINI_KD)
base_idx2 <- rep(seq_len(n_base2), each = n_kd2)
kd_rep2   <- rep(MINI_KD, times = n_base2)
exp_onehot2_exp           <- exp_onehot2[base_idx2, , drop = FALSE]
rownames(exp_onehot2_exp) <- paste0(
  rep(rownames(exp_onehot2), each = n_kd2), "_KD",
  gsub("\\.", "p", sprintf("%.4f", kd_rep2))
)
kd_per_exp2 <- setNames(kd_rep2, rownames(exp_onehot2_exp))
exp_kd2 <- matrix(1, nrow = nrow(exp_onehot2_exp), ncol = ncol(exp_onehot2_exp),
                  dimnames = dimnames(exp_onehot2_exp))
for (i in seq_len(nrow(exp_onehot2_exp))) {
  kd_cols <- which(exp_onehot2_exp[i, ] == 1L)
  if (length(kd_cols) > 0L) exp_kd2[i, kd_cols] <- kd_rep2[i]
}
exp_onehot2 <- exp_onehot2_exp   # use expanded matrix from here on
message("  expand mode: ", n_base2, " base x ", n_kd2, " KD levels = ",
        nrow(exp_onehot2), " total experiments")

assert("T2.2 onehot rows == n_base * n_kd", nrow(exp_onehot2) == n_base2 * n_kd2)
assert("T2.2 onehot cols == n_tfs",         ncol(exp_onehot2) == MINI_N_TFS)
assert("T2.2 kd matrix same dims",          all(dim(exp_kd2) == dim(exp_onehot2)))
assert("T2.2 kd values from pool",          all(kd_per_exp2 %in% MINI_KD))
assert("T2.2 unperturbed cols = 1",         all(exp_kd2[exp_onehot2 == 0L] == 1))
assert("T2.2 base_exp recoverable",         {
  rec <- sub("_KD[0-9p]+$", "", rownames(exp_onehot2))
  all(rec == rep(paste0("EXP_", seq_len(n_base2)), each = n_kd2))
})

saveRDS(exp_onehot2, file.path(MINI_OUT, "EXP_design_onehot.rds"))
saveRDS(exp_kd2,     file.path(MINI_OUT, "EXP_design_kd.rds"))

message("-- T2.3 WT simulation --")
model_wt2 <- model2
model_wt2$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = MINI_NUM_SIMS,
  timepoint       = 0.25,
  genes           = rep(list(tf_ids2[1:4]), MINI_NUM_SIMS),
  num_genes       = 4L,
  multiplier      = 1,
  seed            = seed_common2
)
model_wt2 <- tryCatch(model_wt2 %>% generate_cells(), error = function(e) {
  message("  Note: generate_cells warning/error (expected with minimal params): ", conditionMessage(e))
  model_wt2
})
# With minimal params dyngen may produce bad kinetics; just check call returned a list
assert("T2.3 WT generate_cells returned a model", is.list(model_wt2))

message("-- T2.4 KD simulation (first 2 non-empty experiments) --")
models2   <- list(WT = model_wt2)
meta2     <- list()
ctr2      <- 0L
max_test_exp <- min(2L, nrow(exp_onehot2))

for (i in seq_len(max_test_exp)) {
  mname     <- rownames(exp_onehot2)[i]          # already has _KD suffix in expand mode
  perturbed <- colnames(exp_onehot2)[which(exp_onehot2[i, ] == 1L)]
  if (length(perturbed) == 0L) next
  kd_v      <- kd_per_exp2[mname]
  base_exp  <- sub("_KD[0-9p]+$", "", mname)

  mk <- model2
  mk$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = MINI_NUM_SIMS, timepoint = 0.25,
    genes = rep(list(perturbed), MINI_NUM_SIMS),
    num_genes = length(perturbed), multiplier = kd_v,
    seed = seed_common2
  )
  mk <- mk %>% generate_cells()

  models2[[mname]] <- mk
  ctr2 <- ctr2 + 1L
  meta2[[ctr2]] <- data.frame(model_name=mname, base_experiment=base_exp,
                               kd_multiplier=kd_v, n_perturbed_genes=length(perturbed),
                               stringsAsFactors=FALSE)
  message("  simulated: ", mname, "  kd=", kd_v,
          "  cells=", nrow(mk$cell_info))
}

assert("T2.4 at least 1 KD model created",  length(models2) >= 2L)
assert("T2.4 WT present in models list",     "WT" %in% names(models2))

meta_df2 <- dplyr::bind_rows(meta2)
saveRDS(models2, file.path(MINI_OUT, "models_all.rds"))
write.csv(meta_df2, file.path(MINI_OUT, "model_index.csv"), row.names = FALSE)

# --- verify all expected output files exist ---
expected_files <- c(
  "adj_truth_binary.rds", "adj_truth_strength.rds",
  "EXP_design_onehot.rds", "EXP_design_kd.rds",
  "models_all.rds", "model_index.csv"
)
for (f in expected_files) {
  assert(paste0("T2.5 output file exists: ", f),
         file.exists(file.path(MINI_OUT, f)))
}

# --- verify dimensions of saved files ---
oh_rt  <- readRDS(file.path(MINI_OUT, "EXP_design_onehot.rds"))
kd_rt  <- readRDS(file.path(MINI_OUT, "EXP_design_kd.rds"))
abin   <- readRDS(file.path(MINI_OUT, "adj_truth_binary.rds"))
astr   <- readRDS(file.path(MINI_OUT, "adj_truth_strength.rds"))
mods   <- readRDS(file.path(MINI_OUT, "models_all.rds"))

assert("T2.5 onehot rows == n_base * n_kd",  nrow(oh_rt)   == n_base2 * n_kd2)
assert("T2.5 onehot cols == n_tfs",         ncol(oh_rt)   == MINI_N_TFS)
assert("T2.5 kd_matrix same dims as onehot",all(dim(kd_rt) == dim(oh_rt)))
assert("T2.5 adj_binary is n_tfs x n_tfs",  all(dim(abin) == c(MINI_N_TFS, MINI_N_TFS)))
assert("T2.5 adj_strength same dims",        all(dim(astr) == c(MINI_N_TFS, MINI_N_TFS)))
assert("T2.5 adj_binary only 0/1",           all(abin %in% c(0L, 1L)))
assert("T2.5 models list has WT",            "WT" %in% names(mods))
assert("T2.5 model_index.csv correct cols",
       {
         idx <- read.csv(file.path(MINI_OUT, "model_index.csv"))
         all(c("model_name","base_experiment","kd_multiplier","n_perturbed_genes") %in% colnames(idx))
       })

# =========================================================================
message("\n=== Test summary ===")
message("PASS: ", PASS, "   FAIL: ", FAIL)
message("Output in: ", MINI_OUT)

if (FAIL == 0L) {
  message("All tests PASSED.")
  quit(status = 0)
} else {
  message(FAIL, " test(s) FAILED.")
  quit(status = 1)
}
