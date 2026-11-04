#' Names of the bundled networks
#'
#' @return character vector of network names usable in [load_network()].
#' @export
bundled_networks <- function() {
  list.files(system.file("networks", package = "perturbdyngen"))
}

#' Path of a bundled network, or of a network directory
#'
#' @param network a bundled network name (see [bundled_networks()]) or a directory
#'   written by [build_network()].
#' @return the directory.
#' @export
network_dir <- function(network) {
  if (dir.exists(network)) return(normalizePath(network))
  d <- system.file("networks", network, package = "perturbdyngen")
  if (!nzchar(d)) stop("unknown network '", network, "'; bundled: ",
                       paste(bundled_networks(), collapse = ", "))
  d
}

#' Load a network
#'
#' @param network a bundled network name or a network directory.
#' @return list with the dyngen `model`, the per-cell `seeds`, the random and hub
#'   designs (`design_random`, `design_hub`, and `design_hub_original` when present),
#'   the true network `grn_truth`, and `dir`.
#' @export
load_network <- function(network) {
  d  <- network_dir(network)
  rd <- function(f) { p <- file.path(d, f); if (file.exists(p)) readRDS(p) else NULL }
  list(model = rd("model.rds"), seeds = rd("seeds.rds"), design_random = rd("design_random.rds"),
       design_hub = rd("design_hub.rds"), design_hub_original = rd("design_hub_original.rds"),
       grn_truth = rd("grn_truth.rds"), dir = d)
}

#' Transcription factor ids of a model
#'
#' Experiment designs index transcription factors only.
#'
#' @param model a dyngen model.
#' @export
model_tf_ids <- function(model) {
  fi <- model$feature_info
  if ("is_tf" %in% colnames(fi)) fi$feature_id[fi$is_tf] else fi$feature_id
}

#' Drop the parts of a dyngen model that simulation does not use
#'
#' Removes stored simulations and the gold standard, which make up almost all of a
#' model's size.
#'
#' @param model a dyngen model.
#' @export
slim_model <- function(model) {
  model[setdiff(names(model), c("simulations", "gold_standard", "experiment"))]
}

#' Rescale protein half-lives of a model
#'
#' Updates half-lives, decay rates and the decay-rate parameters of the simulation
#' system, then recomputes `max_protein`.
#'
#' @param model a dyngen model.
#' @param halflife named vector of new half-lives per feature id.
#' @export
update_protein_halflife <- function(model, halflife) {
  fi <- model$feature_info
  for (feature_id in names(halflife)) {
    decay <- log(2) / halflife[[feature_id]]
    fi[fi$feature_id == feature_id, "protein_halflife"] <- halflife[[feature_id]]
    fi[fi$feature_id == feature_id, "protein_decay_rate"] <- decay
    model$simulation_system$parameters[paste0("protein_decay_rate_", feature_id)] <- as.numeric(decay)
  }
  fi$max_protein <- with(fi, translation_rate / protein_decay_rate * max_mrna)
  model$feature_info <- fi
  model
}

#' Default settings of a network build
#'
#' @return named list; see the bundled files in `inst/configs/networks` for the meaning
#'   of each field.
#' @export
default_network_config <- function() {
  list(type = "network", template = "legacy", topology = "sampled", topology_file = NULL,
       sparsity = NULL, n_edges = NULL, sdlog = 1.2, frac_zero_out = 0.25, max_out = 100L,
       frac_self_loops = 0.0901, max_repressors_per_gene = 3L, effect_probs = c(0.75, 0.25),
       strength_values = 1, hill_values = c(2, 4), spine_len = 5L, burn_frac = 0.10,
       protein_halflife_scale = 3, n_cells = 100L, n_targets = 0L, hub_top_k = 100L,
       hub_min_upstream = 2L, max_hub_experiments = 19L, max_chunks_per_hub = 2L,
       include_repressors = FALSE, max_nodes = 8L, seed = 2025L)
}

#' Build a network
#'
#' Samples a topology over the template network's genes (or reads an external one),
#' adds burn-in scaffolding, fixes the edge count, repairs reachability and builds the
#' dyngen model with the template's differentiation path, simulation settings and
#' seeds. The random design is the template's; the hub design is rebuilt for the new
#' topology.
#'
#' @param config path of a network config file or a list (see
#'   [default_network_config()]).
#' @param out_dir output directory.
#' @param cache_dir dyngen download cache.
#' @return the output directory, invisibly. It holds `model.rds`, `seeds.rds`,
#'   `design_random.rds`, `design_hub.rds`, `grn_truth.rds`, `module_network.tsv`,
#'   `summary.csv` and, for external topologies, `gene_to_node.csv`.
#' @export
build_network <- function(config, out_dir, cache_dir = file.path(tempdir(), "dyngen_cache")) {
  cfg  <- utils::modifyList(default_network_config(), read_config(config))
  tmpl <- load_network(cfg$template)
  base <- tmpl$model
  bb    <- base$backbone
  nodes <- bb$module_info$module_id
  n     <- length(nodes)
  sp    <- base$simulation_params
  seed  <- as.integer(cfg$seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (cfg$topology == "external") {
    topo <- read_external_topology(topology_path(cfg$topology_file), nodes,
                                   cfg$max_repressors_per_gene, seed = seed)
    n_edges <- nrow(topo$module_network)
    utils::write.csv(topo$gene_map, file.path(out_dir, "gene_to_node.csv"), row.names = FALSE)
  } else {
    n_edges <- if (!is.null(cfg$n_edges)) as.integer(cfg$n_edges) else
      as.integer(round((1 - cfg$sparsity) * n * n))
    topo <- sample_topology(nodes, n_edges = n_edges, sdlog = cfg$sdlog,
                            frac_zero_out = cfg$frac_zero_out, frac_self_loops = cfg$frac_self_loops,
                            max_repressors_per_gene = cfg$max_repressors_per_gene,
                            max_out = as.integer(cfg$max_out), effect_probs = cfg$effect_probs,
                            strength_values = cfg$strength_values, hill_values = cfg$hill_values,
                            seed = seed)
  }
  burn <- build_burn_structure(topo$module_network, nodes, spine_len = cfg$spine_len,
                               burn_frac = cfg$burn_frac,
                               n_self_loops = length(topo$self_loop_nodes), seed = seed)
  spine_edges <- bind_rows(lapply(burn$spines, function(s)
    tibble(from = utils::head(s, -1), to = utils::tail(s, -1))))
  mn <- enforce_edge_budget(burn$module_network, nodes, n_edges, roots = burn$roots,
                            spine_edges = spine_edges, seed = seed)
  mn <- repair_connectivity(mn, nodes, roots = burn$roots, seed = seed)$module_network
  mn <- cap_repressors(mn, cfg$max_repressors_per_gene, seed = seed)

  new_bb <- dyngen::backbone(module_info = burn$module_info, module_network = mn,
                             expression_patterns = bb$expression_patterns)
  new_bb$expression_patterns <- bb$expression_patterns

  init <- dyngen::initialise_model(
    backbone = new_bb, num_cells = as.integer(cfg$n_cells), num_targets = as.integer(cfg$n_targets),
    num_hks = 0L, experiment = base$experiment_params, verbose = FALSE,
    download_cache_dir = cache_dir,
    simulation_params = dyngen::simulation_default(
      census_interval   = sp$census_interval,
      ssa_algorithm     = GillespieSSA2::ssa_etl(tau = sp$ssa_algorithm$params$tau),
      experiment_params = dyngen::simulation_type_wild_type(
        num_simulations = nrow(sp$experiment_params)),
      total_time = sp$total_time))
  init$simulation_params$experiment_params <- sp$experiment_params
  init$simulation_params$burn_time         <- sp$burn_time

  model <- init |>
    dyngen::generate_tf_network() |>
    dyngen::generate_feature_network() |>
    dyngen::generate_kinetics() |>
    dyngen::generate_gold_standard()

  hl <- suppressWarnings(as.numeric(cfg$protein_halflife_scale))
  if (length(hl) == 1 && !is.na(hl) && hl > 0 && hl != 1) {
    model <- update_protein_halflife(
      model, stats::setNames(model$feature_info$protein_halflife * hl, model$feature_info$feature_id))
    attr(model, "protein_halflife_scale") <- hl
  }

  truth  <- network_truth(model)
  tf_ids <- model_tf_ids(model)
  random <- tmpl$design_random
  if (ncol(random) != length(tf_ids)) stop("template random design does not match the network size")
  colnames(random) <- tf_ids
  hub <- hub_design(truth[tf_ids, tf_ids, drop = FALSE], max_nodes = cfg$max_nodes,
                    top_k = cfg$hub_top_k, min_upstream = cfg$hub_min_upstream,
                    include_repressors = cfg$include_repressors,
                    max_experiments = cfg$max_hub_experiments,
                    max_chunks_per_hub = cfg$max_chunks_per_hub)

  saveRDS(slim_model(model), file.path(out_dir, "model.rds"))
  saveRDS(init$simulation_params$experiment_params$seed, file.path(out_dir, "seeds.rds"))
  saveRDS(random, file.path(out_dir, "design_random.rds"))
  saveRDS(hub, file.path(out_dir, "design_hub.rds"))
  saveRDS(truth, file.path(out_dir, "grn_truth.rds"))
  utils::write.table(as.data.frame(mn), file.path(out_dir, "module_network.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  stats_row <- cbind(topology_stats(mn, nodes, burn$burn_ids),
                     n_roots = length(burn$roots), coverage = burn$coverage,
                     n_random_experiments = nrow(random), n_hub_experiments = nrow(hub))
  utils::write.csv(stats_row, file.path(out_dir, "summary.csv"), row.names = FALSE)
  message(sprintf("network: %d edges, sparsity %.4f, %d random + %d hub experiments -> %s",
                  stats_row$n_edges, stats_row$sparsity, nrow(random), nrow(hub), out_dir))
  invisible(out_dir)
}

#' Path of a bundled external topology, or of a topology file
#'
#' @param x a bundled topology name (file stem in `inst/topologies`) or a file path.
#' @export
topology_path <- function(x) {
  if (file.exists(x)) return(normalizePath(x))
  p <- system.file("topologies", paste0(x, ".tsv"), package = "perturbdyngen")
  if (!nzchar(p)) stop("unknown topology '", x, "'")
  p
}
