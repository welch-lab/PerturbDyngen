#' Default settings of a simulated screen
#'
#' @return named list; see the bundled files in `inst/configs/datasets` for the meaning
#'   of each field.
#' @export
default_dataset_config <- function() {
  list(type = "dataset", network = "legacy", hub_design = "pure", perturbation = "knockdown",
       levels = c(0.005, 0.2, 0.4, 0.6, 0.8), experiments = "all", n_sims = 300L, seed = 2025L,
       pre_kd_window = 0.05, snap_time = "last", guide_sigma_s2 = 0.1, guide_nb_theta = 10,
       sigmoid_k = 0.01484, sigmoid_u0 = 143.4, umi_capture_rate = 0.9, act_guide_kd_equiv = 0.2,
       p_rna = 0.05, p_premrna = 0.01, rna_nb_theta = 3)
}

#' Experiment design of a screen
#'
#' Stacks the random and hub designs and repeats every experiment at each perturbation
#' level. Rows are named `EXP_<i>_KD<level>` / `HUB_EXP_<i>_ACT<level>`.
#'
#' @param net a network from [load_network()].
#' @param cfg a dataset config (see [default_dataset_config()]).
#' @return list with the 0/1 `onehot` design, the multiplier design `multiplier`
#'   (1 for unperturbed genes) and the `level` of each experiment.
#' @export
screen_design <- function(net, cfg) {
  tf_ids <- model_tf_ids(net$model)
  name_rows <- function(m, prefix) {
    if (ncol(m) != length(tf_ids)) stop("design has ", ncol(m), " columns but the model has ",
                                        length(tf_ids), " transcription factors")
    m <- unclass(m); attr(m, "hub_of_experiment") <- NULL
    dimnames(m) <- list(if (nrow(m)) paste0(prefix, seq_len(nrow(m))) else NULL, tf_ids)
    m
  }
  random <- name_rows(net$design_random, "EXP_")
  hub <- switch(cfg$hub_design, pure = net$design_hub, original = net$design_hub_original,
                stop("hub_design must be 'pure' or 'original'"))
  if (is.null(hub)) stop("network has no ", cfg$hub_design, " hub design")
  hub <- name_rows(hub, "HUB_EXP_")
  base <- switch(cfg$experiments, random = random, hub = hub, all = rbind(random, hub),
                 stop("experiments must be 'random', 'hub' or 'all'"))

  tag    <- if (cfg$perturbation == "knockdown") "KD" else "ACT"
  levels <- as.numeric(cfg$levels)
  idx    <- rep(seq_len(nrow(base)), each = length(levels))
  lv     <- rep(levels, times = nrow(base))
  onehot <- base[idx, , drop = FALSE]
  rownames(onehot) <- paste0(rownames(base)[idx], "_", tag, gsub("\\.", "p", sprintf("%.4f", lv)))
  mult <- onehot
  for (i in seq_len(nrow(onehot))) {
    on <- onehot[i, ] == 1L
    mult[i, on] <- lv[i]
    mult[i, !on] <- 1
  }
  list(onehot = onehot, multiplier = mult, level = stats::setNames(lv, rownames(onehot)))
}

#' Simulate a Perturb-seq screen
#'
#' For every experiment, simulates `n_sims` cells as wild-type / perturbed pairs
#' ([simulate_experiment()]), then guide UMIs ([simulate_guide_umis()]) and RNA UMIs
#' at the readout snapshot ([simulate_rna_umis()]) for both copies with shared size
#' factors. Each stage is saved per experiment and skipped if its file exists, so an
#' interrupted run resumes. Seeds depend only on the run seed and the experiment's row
#' in the full design, so running a subset of experiments gives the same files.
#'
#' @param config path, bundled config name or list (see [default_dataset_config()]).
#' @param out_dir output directory.
#' @param n_cores experiments simulated in parallel (forked).
#' @param experiments optional experiment names to run; NULL = all.
#' @return the output directory, invisibly. It holds `checkpoints/<exp>.rds` (true
#'   molecule counts), `perturb_seq/individual/<exp>.rds` (guide UMIs),
#'   `perturb_seq/rna_individual/<exp>.rds` (RNA UMIs), `design/`, `model_index.csv`
#'   and `run_manifest.csv`.
#' @export
simulate_screen <- function(config, out_dir, n_cores = 1L, experiments = NULL) {
  cfg <- utils::modifyList(default_dataset_config(), read_config(config))
  cfg$perturbation <- match.arg(cfg$perturbation, c("knockdown", "activation"))
  levels <- as.numeric(cfg$levels)
  if (cfg$perturbation == "knockdown" && any(levels <= 0 | levels >= 1))
    stop("knockdown levels must be in (0, 1)")
  if (cfg$perturbation == "activation" && any(levels <= 1)) stop("activation levels must be > 1")

  net    <- load_network(cfg$network)
  design <- screen_design(net, cfg)
  n_sims <- as.integer(cfg$n_sims)
  if (length(net$seeds) < n_sims) stop("network has ", length(net$seeds), " seeds; n_sims = ", n_sims)
  seeds  <- net$seeds[seq_len(n_sims)]
  sys    <- simulation_setup(net$model, cfg$pre_kd_window)
  seed0  <- as.integer(cfg$seed)
  sigma_s <- sqrt(cfg$guide_sigma_s2)
  act_umi <- if (cfg$perturbation == "activation")
    sigmoid_inverse(cfg$act_guide_kd_equiv, cfg$sigmoid_k, cfg$sigmoid_u0) else NULL

  d <- list(ckpt = file.path(out_dir, "checkpoints"), design = file.path(out_dir, "design"),
            guide = file.path(out_dir, "perturb_seq", "individual"),
            rna = file.path(out_dir, "perturb_seq", "rna_individual"))
  for (p in d) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  saveRDS(design$onehot, file.path(d$design, "EXP_design_onehot.rds"))
  saveRDS(design$multiplier, file.path(d$design, "EXP_design_kd.rds"))
  saveRDS(as.data.frame(net$model$feature_info), file.path(d$design, "feature_info.rds"))
  saveRDS(net$grn_truth, file.path(d$design, "grn_truth.rds"))

  all_exps <- rownames(design$onehot)
  exps <- if (is.null(experiments)) all_exps else intersect(all_exps, experiments)
  if (!length(exps)) stop("no matching experiments")
  message(sprintf("%s: %d experiments x %d cells, %d cores", cfg$network, length(exps), n_sims, n_cores))

  run_one <- function(e) {
    f_ckpt  <- file.path(d$ckpt, paste0(e, ".rds"))
    f_guide <- file.path(d$guide, paste0(e, ".rds"))
    f_rna   <- file.path(d$rna, paste0(e, ".rds"))
    eseed   <- experiment_seed(seed0, match(e, all_exps))
    level   <- design$level[[e]]

    if (!file.exists(f_ckpt)) {
      genes <- colnames(design$onehot)[design$onehot[e, ] == 1L]
      if (!length(genes)) return(invisible(NULL))
      t0  <- proc.time()[["elapsed"]]
      res <- simulate_experiment(genes, level, sys, seeds)
      if (is.null(res)) { message("  [fail] ", e); return(invisible(NULL)) }
      saveRDS(c(list(experiment = e, perturbed_genes = genes), res), f_ckpt)
      message(sprintf("  [done] %s  %.1f min", e, (proc.time()[["elapsed"]] - t0) / 60))
    }

    if (!file.exists(f_guide)) {
      ck <- readRDS(f_ckpt)
      g  <- simulate_guide_umis(n_cells = ck$n_sims, genes = ck$perturbed_genes, kd_level = level,
                                sigma_s = sigma_s, nb_theta = cfg$guide_nb_theta,
                                capture_rate = cfg$umi_capture_rate, seed = eseed,
                                u_target = act_umi, sigmoid_k = cfg$sigmoid_k,
                                sigmoid_u0 = cfg$sigmoid_u0)
      saveRDS(list(model_name = e, perturbed_genes = ck$perturbed_genes, kd_level = level,
                   n_cells = ck$n_sims, u_target = g$u_target, u_raw = g$u_raw,
                   true_umis = g$true_umis, obs_umis = g$obs_umis,
                   kd_continuous = g$kd_continuous, s_cells = g$s_cells), f_guide)
    }

    if (!file.exists(f_rna)) {
      ck <- readRDS(f_ckpt)
      s_cells <- readRDS(f_guide)$s_cells
      si <- readout_index(ck$snap_times, cfg$snap_time)
      set.seed(eseed + 500L)
      kd <- simulate_rna_umis(ck$mat_kd, ck$mol_ids, si, cfg$p_rna, cfg$p_premrna,
                              cfg$rna_nb_theta, s_cells)
      wt <- simulate_rna_umis(ck$mat_wt, ck$mol_ids, si, cfg$p_rna, cfg$p_premrna,
                              cfg$rna_nb_theta, s_cells)
      saveRDS(list(model_name = e, snap_time = ck$snap_times[si], snap_idx = si,
                   obs_spliced_kd = kd$obs_spliced, obs_unspliced_kd = kd$obs_unspliced,
                   obs_spliced_wt = wt$obs_spliced, obs_unspliced_wt = wt$obs_unspliced), f_rna)
    }
    invisible(NULL)
  }

  parallel::mclapply(exps, run_one, mc.cores = n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  write_run_tables(out_dir, cfg, design, exps, act_umi)
  invisible(out_dir)
}

#' Index of the readout snapshot
#'
#' @param snap_times snapshot times.
#' @param snap_time "last" or a time; the nearest snapshot is used.
#' @export
readout_index <- function(snap_times, snap_time = "last") {
  if (identical(as.character(snap_time), "last")) return(length(snap_times))
  t <- suppressWarnings(as.numeric(snap_time))
  if (is.na(t)) stop("snap_time must be 'last' or a number")
  which.min(abs(snap_times - t))
}

#' Write the parameter and index tables of a run
#'
#' @param out_dir run directory.
#' @param cfg dataset config.
#' @param design design from [screen_design()].
#' @param exps experiments of this call.
#' @param act_umi guide UMI target of activation screens, or NULL.
#' @keywords internal
write_run_tables <- function(out_dir, cfg, design, exps, act_umi) {
  ps <- function(...) file.path(out_dir, "perturb_seq", ...)
  all_exps <- rownames(design$onehot)
  has <- function(sub) all_exps[file.exists(file.path(out_dir, sub, paste0(all_exps, ".rds")))]
  done_crn <- has("checkpoints")
  done_gd  <- has(file.path("perturb_seq", "individual"))
  done_rna <- has(file.path("perturb_seq", "rna_individual"))
  umi_target <- function(lv) if (is.null(act_umi)) sigmoid_inverse(lv, cfg$sigmoid_k, cfg$sigmoid_u0) else
    rep(act_umi, length(lv))

  levels <- as.numeric(cfg$levels)
  utils::write.csv(data.frame(kd_level = levels, u_target = umi_target(levels),
                              u_raw = umi_target(levels) / exp(cfg$guide_sigma_s2 / 2)),
                   ps("kd_level_umi_targets.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    parameter = c("kd_levels", "guide_sigma_s2", "guide_nb_theta", "sigmoid_k", "sigmoid_u0",
                  "umi_capture_rate"),
    value = c(paste(levels, collapse = ","), cfg$guide_sigma_s2, cfg$guide_nb_theta, cfg$sigmoid_k,
              cfg$sigmoid_u0, cfg$umi_capture_rate)), ps("guide_params.csv"), row.names = FALSE)

  snap <- if (length(done_rna)) readRDS(ps("rna_individual", paste0(done_rna[1], ".rds")))$snap_time else NA
  utils::write.csv(data.frame(
    parameter = c("p_rna", "p_premrna", "rna_nb_theta", "snap_time_param", "actual_snap_time"),
    value = c(cfg$p_rna, cfg$p_premrna, cfg$rna_nb_theta, cfg$snap_time, snap)),
    ps("rna_params.csv"), row.names = FALSE)

  if (length(done_crn)) {
    lv <- design$level[done_crn]
    utils::write.csv(data.frame(
      model_name = done_crn, base_experiment = sub("_(KD|ACT)[0-9p]+$", "", done_crn),
      kd_level = lv, umi_target = umi_target(lv),
      n_perturbed_genes = rowSums(design$onehot[done_crn, , drop = FALSE]),
      has_guide_umis = done_crn %in% done_gd, has_rna_umis = done_crn %in% done_rna),
      file.path(out_dir, "model_index.csv"), row.names = FALSE)
  }

  flat <- lapply(cfg, function(v) if (length(v) > 1) paste(v, collapse = ",") else
    if (is.null(v)) NA else v)
  man <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), flat,
                    perturbation_tag = if (cfg$perturbation == "knockdown") "kd" else "act",
                    n_total_experiments = length(all_exps), n_crn_complete = length(done_crn),
                    n_guide_complete = length(done_gd), n_rna_complete = length(done_rna),
                    package_version = as.character(utils::packageVersion("perturbdyngen")),
                    stringsAsFactors = FALSE)
  utils::write.csv(man, file.path(out_dir, "run_manifest.csv"), row.names = FALSE)
}
