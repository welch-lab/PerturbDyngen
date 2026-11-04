#' Guide UMI count expected at a knockdown level
#'
#' Inverse of the sigmoid that maps a guide UMI count to knockdown efficiency.
#'
#' @param kd remaining expression fraction, in (0, 1).
#' @param k,u0 sigmoid steepness and midpoint (UMI count at kd = 0.5).
#' @export
sigmoid_inverse <- function(kd, k = 0.01484, u0 = 143.4) u0 + log(1 / kd - 1) / k

#' Knockdown level implied by a guide UMI count
#'
#' @param umi guide UMI count.
#' @param k,u0 sigmoid steepness and midpoint.
#' @export
sigmoid_kd <- function(umi, k = 0.01484, u0 = 143.4) 1 / (1 + exp(k * (umi - u0)))

#' Simulate guide UMI counts for one experiment
#'
#' Each cell gets a library size factor `s ~ LogNormal(0, sigma_s^2)`. Guide UMIs are
#' `NB(mu = s * u_target / E[s], size = nb_theta)` per targeted gene, and observed UMIs
#' are a binomial capture of those. The knockdown level read back from a cell's mean
#' guide UMI is returned as `kd_continuous`.
#'
#' @param n_cells number of cells.
#' @param genes targeted genes.
#' @param kd_level knockdown level that sets the UMI target; ignored if `u_target` is
#'   given.
#' @param sigma_s standard deviation of the log size factor.
#' @param nb_theta negative binomial size.
#' @param capture_rate binomial capture probability.
#' @param seed seed set before drawing, or NULL.
#' @param u_target fixed UMI target (activation screens, where the knockdown sigmoid
#'   does not apply); `kd_continuous` is then NA.
#' @param sigmoid_k,sigmoid_u0 sigmoid parameters.
#' @return list with `true_umis`, `obs_umis` (cells x genes), `kd_continuous`, the size
#'   factors `s_cells`, `u_raw` and `u_target`.
#' @export
simulate_guide_umis <- function(n_cells, genes, kd_level, sigma_s = sqrt(0.1), nb_theta = 10,
                                capture_rate = 0.9, seed = NULL, u_target = NULL,
                                sigmoid_k = 0.01484, sigmoid_u0 = 143.4) {
  if (!is.null(seed)) set.seed(seed)
  fixed_target <- !is.null(u_target)
  if (!fixed_target) {
    if (kd_level <= 0 || kd_level >= 1) stop("kd_level must be in (0, 1)")
    u_target <- sigmoid_inverse(kd_level, sigmoid_k, sigmoid_u0)
  }
  if (u_target <= 0) stop("non-positive guide UMI target; adjust the sigmoid parameters")

  n_genes <- length(genes)
  u_raw   <- u_target / exp(sigma_s^2 / 2)
  cells   <- paste0("cell_", seq_len(n_cells))
  s_cells <- stats::rlnorm(n_cells, meanlog = 0, sdlog = sigma_s)

  true_umis <- matrix(0L, n_cells, n_genes, dimnames = list(cells, genes))
  for (c_idx in seq_len(n_cells))
    true_umis[c_idx, ] <- as.integer(stats::rnbinom(n_genes, size = nb_theta,
                                                    mu = s_cells[c_idx] * u_raw))
  obs_umis <- matrix(as.integer(stats::rbinom(n_cells * n_genes, size = as.integer(true_umis),
                                              prob = capture_rate)),
                     n_cells, n_genes, dimnames = dimnames(true_umis))

  mean_umi <- if (n_genes > 1L) rowMeans(true_umis) else true_umis[, 1L]
  kd_continuous <- if (fixed_target) stats::setNames(rep(NA_real_, n_cells), cells) else
    sigmoid_kd(mean_umi, sigmoid_k, sigmoid_u0)

  list(true_umis = true_umis, obs_umis = obs_umis, kd_continuous = kd_continuous,
       s_cells = s_cells, u_raw = u_raw, u_target = u_target)
}

#' Simulate observed RNA UMI counts at one snapshot
#'
#' For cell c and gene g, `mu = true_count * p * s_c` and the observed count is
#' `min(NB(mu, size = nb_theta), true_count)`. Spliced counts are drawn first, then
#' unspliced counts. `s_cells` are the guide-layer size factors, so guide and RNA
#' library sizes are shared.
#'
#' @param mat molecule x snapshot x cell array of true counts.
#' @param mol_ids molecule ids (rows of `mat`).
#' @param snap_idx snapshot to read out.
#' @param p_rna,p_premrna capture probabilities of spliced and unspliced RNA.
#' @param nb_theta negative binomial size.
#' @param s_cells size factor per cell.
#' @return list with integer cells x genes matrices `obs_spliced` and `obs_unspliced`.
#' @export
simulate_rna_umis <- function(mat, mol_ids, snap_idx, p_rna, p_premrna, nb_theta, s_cells) {
  mrna_idx    <- grepl("^mol_mrna_", mol_ids)
  premrna_idx <- grepl("^mol_premrna_", mol_ids)
  genes   <- sub("^mol_mrna_", "", mol_ids[mrna_idx])
  n_cells <- dim(mat)[3]
  n_genes <- sum(mrna_idx)
  dn <- list(paste0("cell_", seq_len(n_cells)), genes)

  true_s <- mat[mrna_idx, snap_idx, ]
  true_u <- mat[premrna_idx, snap_idx, ]
  mu_s <- t(true_s) * (s_cells * p_rna)
  mu_u <- t(true_u) * (s_cells * p_premrna)

  obs_s <- matrix(pmin(as.integer(stats::rnbinom(n_cells * n_genes, size = nb_theta,
                                                 mu = as.vector(mu_s))),
                       as.integer(t(true_s))),
                  n_cells, n_genes, dimnames = dn)
  obs_u <- matrix(pmin(as.integer(stats::rnbinom(n_cells * n_genes, size = nb_theta,
                                                 mu = as.vector(mu_u))),
                       as.integer(t(true_u))),
                  n_cells, n_genes, dimnames = dn)
  list(obs_spliced = obs_s, obs_unspliced = obs_u)
}

#' Seed of one experiment
#'
#' @param seed run seed.
#' @param exp_index row of the experiment in the full design (1-based).
#' @export
experiment_seed <- function(seed, exp_index) as.integer(seed) + as.integer(exp_index) * 1000L

#' Seed of the RNA noise draw of one experiment at one snapshot
#'
#' The readout snapshot uses `exp_seed + 500`, the draw made during simulation; other
#' snapshots use `exp_seed + 600 + snap_idx`.
#'
#' @param exp_seed experiment seed from [experiment_seed()].
#' @param snap_idx snapshot index.
#' @param readout_idx index of the readout snapshot.
#' @export
rna_noise_seed <- function(exp_seed, snap_idx, readout_idx) {
  if (snap_idx > 399L) stop("more than 399 snapshots would overlap seed blocks")
  if (snap_idx == readout_idx) exp_seed + 500L else exp_seed + 600L + as.integer(snap_idx)
}
