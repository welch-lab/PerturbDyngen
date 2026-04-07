# =============================================================
# simulation_config_claude.R
# Default configuration for Simulation_generalized_claude.R
#
# Usage:
#   Rscript Simulation_generalized_claude.R --config=simulation_config_claude.R
#   Any --key=value CLI argument overrides the matching field here.
#
# Copy and rename this file for each new experimental condition.
# =============================================================

cfg <- list(

  # ---- Gene / cell counts -----------------------------------------
  n_tfs        = 500L,   # backbone TF gene count (backbone nodes)
  n_targets    = 0L,     # extra non-TF target genes; 0 = TF-only
  n_cells      = 100L,   # cells per simulated condition
  num_sims     = 300L,   # simulations per model

  # ---- Backbone graph ---------------------------------------------
  sparsity         = 70,     # as integer percent (70 = 70% sparse)
  force_cycle      = TRUE,   # ensure at least one cycle exists
  force_connected  = FALSE,  # ensure weak connectivity
  deadend_prob     = 0.1,    # prob of building cycle from a dead-end node
  cycle_length_min = 3L,     # min cycle length when forcing a cycle
  cycle_length_max = 5L,     # max cycle length when forcing a cycle
  # max_in / max_out: NULL → 0.9 * n_tfs (very permissive, let sparsity control)
  max_in           = NULL,
  max_out          = NULL,
  basal_value      = 0.2,
  effect_probs     = c(0.75, 0.25),              # activation vs repression
  strength_values  = c(1, 2, 4, 10, 100),
  strength_probs   = c(0.4, 0.3, 0.15, 0.1, 0.05),

  # ---- Expression patterns ----------------------------------------
  target_states       = 6L,     # branching states in expression trajectory
  # max_reg_genes: NULL → ceiling(0.05 * n_tfs)
  max_reg_genes       = NULL,
  max_terminal_states = 2L,
  reg_sign_probs      = c(0.75, 0.25),   # + vs - regulation signs

  # ---- Experiment design ------------------------------------------
  max_nodes_exp    = 8L,   # max genes per experiment (always 8)
  # n_experiments: 0 → auto-find minimal K via find_minimal_fair_ordered_K()
  n_experiments    = 0L,
  hub_top_k        = 5L,   # top-k hub genes for hub-group design
  hub_min_upstream = 3L,   # min upstream regulators to qualify as hub

  # ---- KD effects -------------------------------------------------
  # kd_mode = "expand" : every base experiment is run at ALL kd_effects levels
  #                       → n_base × length(kd_effects) total models
  # kd_mode = "sample" : each base experiment draws ONE random level
  #                       → n_base total models
  kd_mode    = "expand",
  kd_effects = c(0.05, 0.1, 0.15, 0.2, 0.25),

  # ---- Reproducibility / output -----------------------------------
  seed       = 2025L,
  save_every = 20L,   # checkpoint every N models
  # out_dir: "" → auto-derive from n_tfs, n_targets, sparsity
  out_dir    = ""
)
