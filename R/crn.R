#' Paired wild-type / perturbed tau-leaping with common random numbers
#'
#' Advances a wild-type and a perturbed copy of the same cell from a shared state.
#' Each tau step draws one U(0, 1) per reaction and both copies fire
#' `qpois(u, propensity * tau)` times, so a reaction whose propensity is the same in
#' both copies fires the same number of times. Genes upstream of the perturbation
#' therefore stay identical in the two copies.
#'
#' @param x0 initial molecule counts, shared by both copies.
#' @param params_wt,params_pert named parameter vectors of the two copies.
#' @param cr compiled reactions from [GillespieSSA2::compile_reactions()].
#' @param tau step size.
#' @param start_time,end_time simulation window.
#' @param census_interval time between recorded snapshots.
#' @return list with `times` and the molecule x snapshot matrices `counts_wt` and
#'   `counts_kd` (the perturbed copy).
#' @export
crn_etl_paired <- function(x0, params_wt, params_pert, cr, tau, start_time, end_time,
                           census_interval) {
  nu       <- cr$state_change
  prop_ptr <- cr$function_pointers
  buf_size <- cr$buffer_size
  n_rxn    <- ncol(nu)
  n_mol    <- nrow(nu)
  propensity <- utils::getFromNamespace("test_propensity_cpp", "GillespieSSA2")

  x_wt <- as.numeric(x0)
  x_kd <- as.numeric(x0)
  t    <- start_time
  next_census <- start_time

  times    <- numeric(0)
  snaps_wt <- vector("list", 0)
  snaps_kd <- vector("list", 0)

  while (t <= end_time + 1e-9) {
    if (t >= next_census - 1e-9) {
      times       <- c(times, t)
      snaps_wt    <- c(snaps_wt, list(x_wt))
      snaps_kd    <- c(snaps_kd, list(x_kd))
      next_census <- next_census + census_interval
    }
    if (t >= end_time - 1e-9) break

    p_wt <- pmax(propensity(prop_ptr, params_wt, buf_size, n_rxn, x_wt, t)$propensity, 0)
    p_kd <- pmax(propensity(prop_ptr, params_pert, buf_size, n_rxn, x_kd, t)$propensity, 0)

    u <- stats::runif(n_rxn)
    k_wt <- stats::qpois(u, p_wt * tau)
    k_kd <- stats::qpois(u, p_kd * tau)

    x_wt <- pmax(0, x_wt + as.numeric(nu %*% k_wt))
    x_kd <- pmax(0, x_kd + as.numeric(nu %*% k_kd))

    t <- t + tau
  }

  list(
    times     = times,
    counts_wt = matrix(unlist(snaps_wt), nrow = n_mol),
    counts_kd = matrix(unlist(snaps_kd), nrow = n_mol)
  )
}

#' Parameters of a perturbed cell
#'
#' Multiplies the transcription rate of every perturbed gene by `multiplier`
#' (< 1 knockdown, > 1 activation). The gene keeps its own regulation.
#'
#' @param params named wild-type parameter vector.
#' @param genes perturbed gene ids.
#' @param multiplier transcription-rate multiplier.
#' @return the modified parameter vector.
#' @export
perturb_params <- function(params, genes, multiplier) {
  for (g in genes) {
    key <- paste0("transcription_rate_", g)
    if (key %in% names(params)) params[key] <- params[key] * multiplier
  }
  params
}

#' Simulate one cell as a wild-type / perturbed pair
#'
#' Runs the burn-in and the unperturbed phase up to the perturbation time with
#' GillespieSSA2's tau-leaping, from the cell's own seed, then continues both copies
#' with [crn_etl_paired()]. Snapshots in the last `pre_window` time units before the
#' perturbation are kept (identical in both copies).
#'
#' @param seed the cell's seed.
#' @param sys simulation settings from [simulation_setup()].
#' @param params_pert parameters of the perturbed copy.
#' @return list with `times`, `counts_wt` and `counts_kd`.
#' @export
simulate_cell <- function(seed, sys, params_pert) {
  set.seed(seed)
  pre <- GillespieSSA2::ssa(
    initial_state     = sys$initial_state,
    reactions         = sys$cr,
    final_time        = sys$burn_time + sys$pert_time,
    params            = sys$params,
    method            = GillespieSSA2::ssa_etl(tau = sys$tau),
    census_interval   = sys$census_interval,
    stop_on_neg_state = FALSE
  )
  x_at_pert <- pre$state[nrow(pre$state), ]

  pre_times <- seq(0, sys$burn_time + sys$pert_time, by = sys$census_interval) - sys$burn_time
  keep      <- which(pre_times >= (sys$pert_time - sys$pre_window) & pre_times < sys$pert_time)
  pre_counts <- t(pre$state[keep, , drop = FALSE])

  out <- crn_etl_paired(
    x0              = x_at_pert,
    params_wt       = sys$params,
    params_pert     = params_pert,
    cr              = sys$cr,
    tau             = sys$tau,
    start_time      = sys$pert_time,
    end_time        = sys$total_time,
    census_interval = sys$census_interval
  )

  list(
    times     = c(pre_times[keep], out$times),
    counts_wt = cbind(pre_counts, out$counts_wt),
    counts_kd = cbind(pre_counts, out$counts_kd)
  )
}

#' Compile a network model for simulation
#'
#' The perturbation is applied at a quarter of the total time, after the burn-in.
#'
#' @param model a dyngen model (see [load_network()]).
#' @param pre_window time before the perturbation that is recorded, as a fraction of
#'   the total time.
#' @return list of compiled reactions, initial state, parameters and timing.
#' @export
simulation_setup <- function(model, pre_window = 0.05) {
  ss <- model$simulation_system
  sp <- model$simulation_params
  total_time <- sp$total_time
  list(
    cr              = GillespieSSA2::compile_reactions(reactions = ss$reactions,
                                                       state_ids = ss$molecule_ids,
                                                       params    = ss$parameters),
    initial_state   = ss$initial_state,
    params          = ss$parameters,
    molecule_ids    = ss$molecule_ids,
    burn_time       = sp$burn_time,
    total_time      = total_time,
    tau             = sp$ssa_algorithm$params$tau,
    census_interval = sp$census_interval,
    pert_time       = 0.25 * total_time,
    pre_window      = pre_window * total_time
  )
}

#' Simulate all cells of one experiment
#'
#' @param genes perturbed gene ids.
#' @param multiplier transcription-rate multiplier.
#' @param sys simulation settings from [simulation_setup()].
#' @param seeds per-cell seeds; cell i uses `seeds[i]`, so the first cells of a run
#'   do not depend on how many cells are simulated.
#' @return list with `snap_times`, `mol_ids`, `n_sims` and the molecule x snapshot x
#'   cell arrays `mat_wt` and `mat_kd`, or NULL if every cell failed.
#' @export
simulate_experiment <- function(genes, multiplier, sys, seeds) {
  params_pert <- perturb_params(sys$params, genes, multiplier)
  cells <- lapply(seq_along(seeds), function(s) {
    tryCatch(simulate_cell(seeds[s], sys, params_pert),
             error = function(e) { message("    cell ", s, " failed: ", conditionMessage(e)); NULL })
  })
  cells <- Filter(Negate(is.null), cells)
  if (!length(cells)) return(NULL)

  n_mol  <- nrow(cells[[1]]$counts_wt)
  n_snap <- length(cells[[1]]$times)
  mat_wt <- array(NA_real_, c(n_mol, n_snap, length(cells)))
  mat_kd <- array(NA_real_, c(n_mol, n_snap, length(cells)))
  for (s in seq_along(cells)) {
    mat_wt[, , s] <- cells[[s]]$counts_wt
    mat_kd[, , s] <- cells[[s]]$counts_kd
  }
  list(snap_times = cells[[1]]$times, mol_ids = sys$molecule_ids, n_sims = length(cells),
       mat_wt = mat_wt, mat_kd = mat_kd)
}
