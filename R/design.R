#' Random experiment design in which every ordered gene pair is separated
#'
#' Assigns each of `n` genes a distinct set of `l` experiments out of `K`, with
#' experiment usage balanced, where `l` is the largest value with
#' `choose(K, l - 1) < n <= choose(K, l)`. Then for most gene pairs (i, j), some
#' experiment perturbs i without j and another perturbs j without i.
#'
#' @param n number of genes.
#' @param K number of experiments.
#' @param max_nodes_per_experiment maximum genes per experiment.
#' @param seed seed.
#' @return list with the K x n 0/1 `experiment_matrix`, `index_sets`, `l`,
#'   `uncovered_pairs` and `covered_percent`.
#' @export
fair_ordered <- function(n, K, max_nodes_per_experiment = 8, seed = 123) {
  set.seed(seed)
  valid_l <- vapply(seq_len(K), function(l) choose(K, l - 1) < n && n <= choose(K, l), logical(1))
  if (!any(valid_l)) stop("no l satisfies choose(K, l - 1) < n <= choose(K, l)")
  l <- max(which(valid_l))

  index_sets <- sample(utils::combn(K, l, simplify = FALSE), n)
  for (iter in seq_len(1000)) {
    freqs <- table(factor(unlist(index_sets), levels = seq_len(K)))
    most  <- as.integer(names(freqs)[which.max(freqs)])
    least <- as.integer(names(freqs)[which.min(freqs)])
    if ((freqs[most] - freqs[least]) <= 1) break
    for (A in utils::combn(setdiff(seq_len(K), c(most, least)), l - 1, simplify = FALSE)) {
      new_M <- sort(c(most, A))
      new_m <- sort(c(least, A))
      if (any(sapply(index_sets, function(s) identical(s, new_M))) &&
          !any(sapply(index_sets, function(s) identical(s, new_m)))) {
        index_sets[[which(sapply(index_sets, function(s) identical(s, new_M)))[1]]] <- new_m
        break
      }
    }
  }

  exp_matrix <- matrix(0, nrow = K, ncol = n)
  for (i in seq_len(n)) exp_matrix[index_sets[[i]], i] <- 1
  if (any(rowSums(exp_matrix) > max_nodes_per_experiment))
    stop("experiments exceed max_nodes_per_experiment; increase K")

  all_pairs <- utils::combn(n, 2, simplify = FALSE)
  uncovered <- Filter(function(p) {
    !(any(exp_matrix[, p[1]] == 1 & exp_matrix[, p[2]] == 0) &&
        any(exp_matrix[, p[2]] == 1 & exp_matrix[, p[1]] == 0))
  }, all_pairs)
  list(experiment_matrix = exp_matrix, index_sets = index_sets, l = l,
       uncovered_pairs = uncovered, covered_percent = 1 - length(uncovered) / length(all_pairs))
}

#' True gene regulatory network of a model
#'
#' @param model a dyngen model.
#' @return genes x genes matrix, `truth[from, to] = effect * strength`.
#' @export
network_truth <- function(model) {
  nodes <- model$feature_info$feature_id
  fn    <- model$feature_network
  fn    <- fn[fn$from %in% nodes & fn$to %in% nodes, , drop = FALSE]
  adj   <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (i in seq_len(nrow(fn))) adj[fn$from[i], fn$to[i]] <- fn$effect[i] * fn$strength[i]
  adj
}

#' Hub experiment design
#'
#' Ranks genes by (signed) in-degree and, for each of the `top_k` highest, perturbs its
#' upstream regulators in one experiment, split into chunks of `max_nodes` when there
#' are more. Each experiment targets the regulators of a single hub. Hubs with fewer
#' than `min_upstream` regulators are skipped, and at most `max_experiments`
#' experiments are kept in hub rank order.
#'
#' @param adj regulator x target truth matrix.
#' @param max_nodes maximum genes per experiment.
#' @param top_k number of hubs considered.
#' @param min_upstream minimum regulators per hub (and per chunk).
#' @param include_repressors FALSE uses activating edges only (`adj == 1`); TRUE uses
#'   any edge.
#' @param max_experiments maximum number of experiments; NULL = no limit.
#' @param max_chunks_per_hub maximum experiments per hub; NULL = no limit.
#' @return experiments x genes 0/1 matrix with attribute `hub_of_experiment`.
#' @export
hub_design <- function(adj, max_nodes = 8L, top_k = 100L, min_upstream = 2L,
                       include_repressors = FALSE, max_experiments = 19L,
                       max_chunks_per_hub = 2L) {
  a <- if (include_repressors) (adj != 0) * 1 else adj
  top_k <- min(top_k, ncol(a))
  hubs <- order(colSums(a), decreasing = TRUE)[seq_len(top_k)]
  sets <- list()
  for (h in hubs) {
    up <- which(a[, h] == 1)
    if (length(up) < min_upstream) next
    chunks <- if (length(up) > max_nodes) split(up, ceiling(seq_along(up) / max_nodes)) else list(up)
    chunks <- Filter(function(ch) length(ch) >= min_upstream, chunks)
    if (!is.null(max_chunks_per_hub)) chunks <- utils::head(chunks, max_chunks_per_hub)
    for (ch in chunks) sets[[length(sets) + 1L]] <- list(genes = ch, hub = h)
  }
  if (!is.null(max_experiments)) sets <- utils::head(sets, max_experiments)
  m <- matrix(0L, length(sets), ncol(adj), dimnames = list(NULL, colnames(adj)))
  for (i in seq_along(sets)) m[i, sets[[i]]$genes] <- 1L
  if (length(sets)) rownames(m) <- paste0("HUB_EXP_", seq_along(sets))
  attr(m, "hub_of_experiment") <- colnames(adj)[vapply(sets, `[[`, integer(1), "hub")]
  m
}
