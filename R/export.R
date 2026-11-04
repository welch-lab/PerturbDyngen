#' Bind reticulate to a Python with anndata
#'
#' Uses `options(perturbdyngen.python)` or the `PERTURBDYNGEN_PYTHON` environment
#' variable when set, otherwise reticulate's default Python.
#'
#' @keywords internal
need_anndata <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("export_h5ad needs the R package 'reticulate'")
  py <- getOption("perturbdyngen.python", Sys.getenv("PERTURBDYNGEN_PYTHON", ""))
  if (nzchar(py) && !reticulate::py_available(initialize = FALSE))
    reticulate::use_python(path.expand(py), required = TRUE)
  if (!reticulate::py_module_available("anndata")) stop("the Python used by reticulate has no anndata")
  invisible(TRUE)
}

#' Read the shared tables of a simulated screen
#'
#' @param run_dir output directory of [simulate_screen()].
#' @return list with the designs, model index, manifest, RNA parameters, feature info,
#'   true network and the label of perturbed cells ("KD" or "ACT").
#' @export
read_run <- function(run_dir) {
  f <- function(...) file.path(run_dir, ...)
  if (!file.exists(f("design", "EXP_design_onehot.rds"))) stop("not a simulated screen: ", run_dir)
  csv <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
  rd  <- function(p) if (file.exists(p)) readRDS(p) else NULL
  run <- list(out_dir = run_dir, design_onehot = readRDS(f("design", "EXP_design_onehot.rds")),
              design_kd = rd(f("design", "EXP_design_kd.rds")), model_index = csv(f("model_index.csv")),
              manifest = csv(f("run_manifest.csv")), rna_params = csv(f("perturb_seq", "rna_params.csv")),
              guide_params = csv(f("perturb_seq", "guide_params.csv")),
              feature_info = rd(f("design", "feature_info.rds")), grn_truth = rd(f("design", "grn_truth.rds")))
  run$pert_label <- if (identical(run$manifest$perturbation[1], "activation")) "ACT" else "KD"
  run
}

#' Time-window label of a snapshot
#'
#' Windows are right-closed bins of `width` time units: group g covers (g * width,
#' (g + 1) * width].
#'
#' @param t snapshot time.
#' @param width bin width.
#' @export
timepoint_group <- function(t, width = 20) {
  as.integer(pmax(ceiling(round(t, 6) / width - 1e-9) - 1L, 0L))
}

#' Molecule x snapshot x cell array to a (cell, snapshot) x gene matrix
#'
#' Rows run over snapshots fastest, then cells.
#'
#' @param arr molecule x snapshot x cell array.
#' @param mol_ids molecule ids.
#' @param prefix molecule prefix, `"mol_mrna_"` or `"mol_premrna_"`.
#' @param snap_idx snapshots to keep.
#' @keywords internal
flatten_counts <- function(arr, mol_ids, prefix, snap_idx) {
  keep <- grep(paste0("^", prefix), mol_ids)
  a <- aperm(arr[keep, snap_idx, , drop = FALSE], c(2L, 3L, 1L))
  m <- matrix(as.integer(a), nrow = dim(a)[1] * dim(a)[2], ncol = dim(a)[3])
  colnames(m) <- sub(paste0("^", prefix), "", mol_ids[keep])
  m
}

#' Sparse block with one design row repeated n times
#'
#' @param dvec 0/1 vector over genes.
#' @param n number of rows.
#' @keywords internal
design_block <- function(dvec, n) {
  j <- which(dvec != 0)
  Matrix::sparseMatrix(i = rep(seq_len(n), times = length(j)), j = rep(j, each = n), x = 1,
                       dims = c(n, length(dvec)))
}

#' Sparse guide UMI block over all genes, each cell repeated once per snapshot
#'
#' @param gm cells x targets guide matrix, or NULL for zeros.
#' @param genes all genes.
#' @param each snapshots per cell.
#' @param n_cells number of cells.
#' @keywords internal
guide_block <- function(gm, genes, each, n_cells) {
  n <- n_cells * each
  if (is.null(gm)) return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                               dims = c(n, length(genes))))
  rows <- gm[rep(seq_len(nrow(gm)), each = each), , drop = FALSE]
  Matrix::sparseMatrix(i = rep(seq_len(n), times = ncol(gm)), j = rep(match(colnames(gm), genes), each = n),
                       x = as.numeric(rows), dims = c(n, length(genes)))
}

#' Row-bind count blocks as CSR (sparse) or dense integer
#'
#' @param blocks list of matrices.
#' @param sparse TRUE for CSR.
#' @keywords internal
bind_blocks <- function(blocks, sparse) {
  if (!sparse) {
    m <- do.call(rbind, blocks); storage.mode(m) <- "integer"
  } else {
    m <- methods::as(do.call(rbind, lapply(blocks, function(b) methods::as(b, "CsparseMatrix"))),
                     "RsparseMatrix")
  }
  dimnames(m) <- NULL
  m
}

#' Export a simulated screen to AnnData
#'
#' One row per cell and snapshot, for every experiment and (once) the wild type:
#'
#' * `X`, `layers$spliced`: noisy spliced UMIs (CSR, int32)
#' * `layers$unspliced`: noisy unspliced UMIs (CSR, int32)
#' * `layers$spliced_true`, `layers$unspliced_true`: true molecule counts (int32)
#' * `obsm$perturb_design`: perturbed genes; `obsm$guide_obs_umis`,
#'   `obsm$guide_true_umis`: guide UMIs over all genes
#' * `uns$grn_truth`: the true network, `truth[from, to] = effect * strength`
#'
#' RNA noise is drawn at every snapshot with the cell's size factor; at the readout
#' snapshot the draw is the one made by [simulate_screen()] (checked when `verify`),
#' and other snapshots use [rna_noise_seed()]. `obs$readout` marks the readout rows.
#' Requires reticulate and a Python with anndata.
#'
#' @param run_dir output directory of [simulate_screen()].
#' @param file h5ad file to write.
#' @param tag dataset label for `obs$tag`.
#' @param experiments experiments to include; NULL = all.
#' @param snap_times snapshot times to include; NULL = all.
#' @param wt "single" (one wild-type copy), "paired" (each experiment's own) or "none".
#' @param guides include guide UMIs in `obsm`.
#' @param verify check the regenerated readout noise of the first experiment against
#'   the stored RNA UMIs.
#' @param max_cells refuse to write more rows than this.
#' @return the file path, invisibly.
#' @export
export_h5ad <- function(run_dir, file, tag = basename(run_dir), experiments = NULL, snap_times = NULL,
                        wt = c("single", "paired", "none"), guides = TRUE, verify = TRUE,
                        max_cells = 2.5e6) {
  wt  <- match.arg(wt)
  run <- read_run(run_dir)
  pth <- function(sub, e) file.path(run_dir, sub, paste0(e, ".rds"))
  ck  <- function(e) pth("checkpoints", e)
  gd  <- function(e) pth(file.path("perturb_seq", "individual"), e)

  all_exps <- rownames(run$design_onehot)
  exps <- if (is.null(experiments)) all_exps else experiments
  exps <- exps[exps %in% all_exps & file.exists(ck(exps)) & file.exists(gd(exps))]
  if (!length(exps)) stop("no experiments with both molecule counts and guide UMIs")

  man   <- run$manifest
  p_rna <- as.numeric(man$p_rna); p_pre <- as.numeric(man$p_premrna)
  theta <- as.numeric(man$rna_nb_theta); seed0 <- as.integer(man$seed)

  first <- readRDS(ck(exps[1]))
  raw_t <- first$snap_times
  all_t <- round(raw_t, 6)
  rp    <- run$rna_params
  ro_t  <- suppressWarnings(as.numeric(rp$value[rp$parameter == "actual_snap_time"]))
  readout_idx <- if (length(ro_t) == 1 && is.finite(ro_t)) which.min(abs(raw_t - ro_t)) else length(raw_t)
  s_idx <- if (is.null(snap_times)) seq_along(all_t) else
    sort(unique(vapply(snap_times, function(t) which.min(abs(all_t - t)), integer(1))))
  times <- all_t[s_idx]
  S <- length(s_idx)
  N <- first$n_sims
  pert_time <- 0.25 * max(all_t)

  n_rows <- S * N * (length(exps) + switch(wt, single = 1L, paired = length(exps), none = 0L))
  if (n_rows > max_cells)
    stop(sprintf("%s rows > max_cells = %s; subset experiments or snap_times",
                 format(n_rows, big.mark = ","), format(max_cells, big.mark = ",")))

  genes  <- sub("^mol_mrna_", "", grep("^mol_mrna_", first$mol_ids, value = TRUE))
  G      <- length(genes)
  kd_lv  <- run$model_index$kd_level[match(exps, run$model_index$model_name)]
  design <- run$design_onehot[exps, genes, drop = FALSE]

  NS <- NU <- TS <- TU <- PD <- GO <- GT <- OBS <- list()
  emit <- function(key, cond, e, kd, dvec, ns, nu, arr, mol_ids, s_c, kdc, gob, gtr) {
    n     <- S * N
    sim_i <- rep(seq_len(N), each = S)
    tt    <- rep(times, times = N)
    single_wt <- cond == "WT" && wt == "single"
    NS[[key]] <<- methods::as(matrix(ns, n, G), "CsparseMatrix")
    NU[[key]] <<- methods::as(matrix(nu, n, G), "CsparseMatrix")
    TS[[key]] <<- flatten_counts(arr, mol_ids, "mol_mrna_", s_idx)[, genes, drop = FALSE]
    TU[[key]] <<- flatten_counts(arr, mol_ids, "mol_premrna_", s_idx)[, genes, drop = FALSE]
    PD[[key]] <<- design_block(if (cond == "KD") dvec else 0L * dvec, n)
    if (guides) {
      GO[[key]] <<- guide_block(if (cond == "KD") gob else NULL, genes, S, N)
      GT[[key]] <<- guide_block(if (cond == "KD") gtr else NULL, genes, S, N)
    }
    OBS[[key]] <<- data.frame(
      tag = tag, model_name = if (single_wt) "WT" else e,
      base_experiment = if (single_wt) "WT" else sub("_(KD|ACT)[0-9p]+$", "", e),
      planned_experiment = e, design_type = ifelse(grepl("^HUB_", e), "hub", "random"),
      condition = cond, kd_level = if (cond == "KD") kd else 1, planned_kd_level = kd,
      n_perturbed_genes = if (cond == "KD") sum(dvec) else 0L, simulation_i = sim_i,
      sim_index = sim_i, sim_time = tt, timepoint_group = timepoint_group(tt),
      post_kd = tt >= pert_time, readout = tt == all_t[readout_idx],
      size_factor = rep(as.numeric(s_c), each = S),
      kd_continuous = if (cond == "KD") rep(as.numeric(kdc), each = S) else NA_real_,
      cell_id = paste0("cell_", sim_i, "_t", tt), stringsAsFactors = FALSE)
  }

  for (j in seq_along(exps)) {
    e   <- exps[j]
    r   <- if (j == 1L) first else readRDS(ck(e))
    gr  <- readRDS(gd(e))
    s_c <- gr$s_cells
    eseed <- experiment_seed(seed0, match(e, all_exps))
    ks <- ku <- ws <- wu <- array(0L, c(S, N, G))
    for (k in seq_len(S)) {
      si <- s_idx[k]
      set.seed(rna_noise_seed(eseed, si, readout_idx))
      nk <- simulate_rna_umis(r$mat_kd, r$mol_ids, si, p_rna, p_pre, theta, s_c)
      nw <- simulate_rna_umis(r$mat_wt, r$mol_ids, si, p_rna, p_pre, theta, s_c)
      ks[k, , ] <- nk$obs_spliced[, genes];  ku[k, , ] <- nk$obs_unspliced[, genes]
      ws[k, , ] <- nw$obs_spliced[, genes];  wu[k, , ] <- nw$obs_unspliced[, genes]
      if (verify && j == 1L && si == readout_idx) {
        ref <- readRDS(pth(file.path("perturb_seq", "rna_individual"), e))
        ok  <- identical(unname(nk$obs_spliced), unname(ref$obs_spliced_kd)) &&
          identical(unname(nk$obs_unspliced), unname(ref$obs_unspliced_kd)) &&
          identical(unname(nw$obs_spliced), unname(ref$obs_spliced_wt)) &&
          identical(unname(nw$obs_unspliced), unname(ref$obs_unspliced_wt))
        if (!ok) stop("regenerated readout noise for ", e, " does not match the stored RNA UMIs")
      }
    }
    dvec <- as.integer(design[e, ])
    emit(paste0(e, "__KD"), "KD", e, kd_lv[j], dvec, ks, ku, r$mat_kd, r$mol_ids, s_c,
         gr$kd_continuous, gr$obs_umis, gr$true_umis)
    if (wt == "paired" || (wt == "single" && j == 1L))
      emit(if (wt == "single") "WT" else paste0(e, "__WT"), "WT", e, kd_lv[j], dvec, ws, wu,
           r$mat_wt, r$mol_ids, s_c, NULL, NULL, NULL)
    rm(r, ks, ku, ws, wu)
    if (j %% 10L == 0L) invisible(gc())
  }

  obs <- do.call(rbind, OBS)
  obs$condition[obs$condition == "KD"] <- run$pert_label
  rownames(obs) <- make.unique(paste(obs$model_name, obs$condition, obs$cell_id, sep = "__"))
  obsm <- list(perturb_design = do.call(rbind, PD))
  if (guides) { obsm$guide_obs_umis <- do.call(rbind, GO); obsm$guide_true_umis <- do.call(rbind, GT) }
  layers <- list(spliced = bind_blocks(NS, TRUE), unspliced = bind_blocks(NU, TRUE),
                 spliced_true = bind_blocks(TS, FALSE), unspliced_true = bind_blocks(TU, FALSE))
  rm(NS, NU, TS, TU); invisible(gc())

  fi  <- run$feature_info
  var <- if (is.null(fi)) data.frame(row.names = genes, feature_id = genes) else {
    rownames(fi) <- fi$feature_id
    fi[genes, setdiff(names(fi), c("mol_premrna", "mol_mrna", "mol_protein", "color")), drop = FALSE]
  }
  df <- function(m) if (is.null(m)) NULL else as.data.frame(m, stringsAsFactors = FALSE)
  uns <- list(tag = tag, conversion_mode = "full", perturb_design_columns = genes,
              exp_design_onehot = df(run$design_onehot[, genes, drop = FALSE]),
              exp_design_kd = if (!is.null(run$design_kd)) df(run$design_kd[, genes, drop = FALSE]),
              run_manifest = man, rna_params = run$rna_params, guide_params = run$guide_params,
              model_index = run$model_index,
              grn_truth = if (!is.null(run$grn_truth)) df(run$grn_truth[genes, genes, drop = FALSE]),
              wt_mode = wt, kd_time = pert_time, snap_times = times, readout_time = all_t[readout_idx],
              noise_params = data.frame(p_rna = p_rna, p_premrna = p_pre, rna_nb_theta = theta,
                                        seed = seed0))
  uns <- uns[!vapply(uns, is.null, logical(1))]

  need_anndata()
  ad_py  <- reticulate::import("anndata", convert = FALSE)
  counts <- function(m) {
    if (inherits(m, "sparseMatrix")) m <- methods::as(m, "RsparseMatrix")
    dimnames(m) <- NULL
    reticulate::r_to_py(m)$astype("int32")
  }
  L  <- lapply(layers, counts)
  ad <- ad_py$AnnData(X = L$spliced, obs = reticulate::r_to_py(obs), var = reticulate::r_to_py(var),
                      obsm = reticulate::dict(lapply(obsm, counts)), layers = reticulate::dict(L),
                      uns = reticulate::r_to_py(uns))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ad$write_h5ad(file, compression = "gzip")
  message(sprintf("wrote %s (%d rows x %d genes)", file, nrow(obs), G))
  invisible(file)
}
