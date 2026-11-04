#' Draw a log-normal out-degree sequence with an exact edge budget
#'
#' Degrees of the active regulators are log-normal, scaled to the budget, clipped to
#' `[1, max_out]` and reconciled proportionally so they sum to `n_edges`. A share
#' `frac_zero_out` of nodes regulates nothing; it is raised when the budget is smaller
#' than the number of active nodes.
#'
#' @param n_nodes number of nodes.
#' @param n_edges edge budget.
#' @param sdlog spread on the log scale; higher gives a heavier hub tail.
#' @param frac_zero_out share of nodes with out-degree 0.
#' @param max_out out-degree cap (clips the distribution).
#' @param seed seed, or NULL to use the current RNG state.
#' @return integer out-degrees in random node order.
#' @export
sample_out_degrees <- function(n_nodes, n_edges, sdlog = 1.0, frac_zero_out = 0.0,
                               max_out = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(max_out)) max_out <- n_nodes - 1L
  max_out <- min(max_out, n_nodes - 1L)

  n_zero   <- as.integer(round(frac_zero_out * n_nodes))
  n_active <- n_nodes - n_zero
  if (n_active < 1L) stop("frac_zero_out leaves no active regulators")
  if (n_edges > n_active * max_out)
    stop("edge budget ", n_edges, " exceeds capacity ", n_active * max_out,
         " (raise max_out or lower frac_zero_out)")
  if (n_edges < n_active) {
    warning("edge budget ", n_edges, " < active nodes ", n_active,
            "; raising frac_zero_out from ", round(frac_zero_out, 3), " to ",
            round(1 - n_edges / n_nodes, 3))
    n_active <- n_edges
    n_zero   <- n_nodes - n_active
  }

  raw <- stats::rlnorm(n_active, meanlog = 0, sdlog = sdlog)
  deg <- raw / sum(raw) * n_edges
  deg <- pmin(pmax(round(deg), 1L), max_out)

  for (iter in seq_len(200L)) {
    diff <- n_edges - sum(deg)
    if (diff == 0L) break
    if (diff > 0L) {
      idx <- which(deg < max_out)
      if (!length(idx)) break
      share <- deg[idx] / sum(deg[idx])
      add   <- as.integer(round(share * diff))
      deg[idx] <- pmin(deg[idx] + add, max_out)
    } else {
      idx <- which(deg > 1L)
      if (!length(idx)) break
      share <- deg[idx] / sum(deg[idx])
      sub   <- as.integer(round(share * abs(diff)))
      deg[idx] <- pmax(deg[idx] - sub, 1L)
    }
  }
  diff <- n_edges - sum(deg)
  if (diff > 0L) {
    idx <- rep(which(deg < max_out), length.out = diff)
    for (i in idx) deg[i] <- min(deg[i] + 1L, max_out)
  } else if (diff < 0L) {
    idx <- rep(which(deg > 1L), length.out = abs(diff))
    for (i in idx) deg[i] <- max(deg[i] - 1L, 1L)
  }
  if (sum(deg) != n_edges) warning("degree sequence sums to ", sum(deg), " vs budget ", n_edges)

  sample(c(deg, rep(0L, n_zero)))
}

#' Cap the number of repressive inputs per gene
#'
#' Repressive edges beyond `max_repressors` on one target are chosen at random and
#' made activating.
#'
#' @param module_network edge table with `from`, `to`, `effect`.
#' @param max_repressors maximum repressive in-edges per target; NULL = no cap.
#' @param seed seed, or NULL.
#' @return the edge table.
#' @export
cap_repressors <- function(module_network, max_repressors = 3L, seed = NULL) {
  if (is.null(max_repressors)) return(module_network)
  if (!is.null(seed)) set.seed(seed)
  mn  <- module_network
  idx <- which(mn$effect == -1L)
  if (!length(idx)) return(mn)
  flip <- unlist(lapply(split(idx, mn$to[idx]), function(ix)
    if (length(ix) > max_repressors) sample(ix, length(ix) - max_repressors) else integer(0)),
    use.names = FALSE)
  if (length(flip)) mn$effect[flip] <- 1L
  mn
}

#' Sample a directed regulatory network
#'
#' Out-degrees come from [sample_out_degrees()]; targets are drawn uniformly among
#' nodes below the in-degree cap. Edge signs, strengths and Hill coefficients are
#' drawn independently, repressors are capped with [cap_repressors()], and
#' self-loops (autoregulation) are added to a share `frac_self_loops` of nodes,
#' preferring the highest out-degree nodes.
#'
#' @param node_ids node names.
#' @param n_edges edge budget (self-loops excluded).
#' @param sdlog,frac_zero_out,max_out degree settings, see [sample_out_degrees()].
#' @param frac_self_loops share of nodes with a self-loop.
#' @param max_repressors_per_gene repressor cap, see [cap_repressors()].
#' @param self_loop_effect,self_loop_hill,self_loop_strength self-loop attributes.
#' @param self_loop_prefer_hubs put self-loops on the highest out-degree nodes.
#' @param max_in in-degree cap; NULL = none.
#' @param effect_probs probabilities of activation and repression.
#' @param strength_values,strength_probs edge strengths and their probabilities.
#' @param hill_values,hill_probs Hill coefficients and their probabilities.
#' @param seed seed, or NULL.
#' @return list with the edge table `module_network` (`from`, `to`, `effect`,
#'   `strength`, `hill`) and `self_loop_nodes`.
#' @export
sample_topology <- function(node_ids, n_edges, sdlog = 1.0, frac_zero_out = 0.0,
                            frac_self_loops = 0.0, max_repressors_per_gene = 3L,
                            self_loop_effect = 1L, self_loop_hill = 4, self_loop_strength = 1,
                            self_loop_prefer_hubs = TRUE, max_out = NULL, max_in = NULL,
                            effect_probs = c(0.75, 0.25), strength_values = 1,
                            strength_probs = NULL, hill_values = c(2, 4), hill_probs = NULL,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(node_ids)
  active <- node_ids

  degs <- sample_out_degrees(n_nodes = length(active), n_edges = n_edges, sdlog = sdlog,
                             frac_zero_out = frac_zero_out, max_out = max_out)
  names(degs) <- sample(active)

  target_pool <- active
  cap_in      <- if (is.null(max_in)) n else max_in
  in_count    <- stats::setNames(integer(length(target_pool)), target_pool)

  edges <- list()
  for (s in names(degs)) {
    d <- degs[[s]]
    if (d < 1L) next
    avail <- target_pool[target_pool != s & in_count[target_pool] < cap_in]
    if (!length(avail)) next
    d  <- min(d, length(avail))
    tg <- if (length(avail) > d) sample(avail, d) else avail
    in_count[tg] <- in_count[tg] + 1L
    edges[[length(edges) + 1L]] <- tibble(from = s, to = tg)
  }
  mn <- bind_rows(edges) |> filter(from != to) |> distinct(from, to)

  mn$effect   <- sample(c(1L, -1L), nrow(mn), replace = TRUE, prob = effect_probs)
  mn$strength <- if (length(strength_values) == 1L) strength_values else
    sample(strength_values, nrow(mn), replace = TRUE, prob = strength_probs)
  mn$hill     <- if (length(hill_values) == 1L) hill_values else
    sample(hill_values, nrow(mn), replace = TRUE, prob = hill_probs)

  mn <- cap_repressors(mn, max_repressors_per_gene)

  self_nodes <- character(0)
  n_self <- as.integer(round(frac_self_loops * n))
  if (n_self > 0L && length(active)) {
    if (self_loop_prefer_hubs) {
      od <- as.integer(table(factor(mn$from, levels = active)))
      self_nodes <- utils::head(active[order(od, decreasing = TRUE)], min(n_self, length(active)))
    } else {
      self_nodes <- sample(active, min(n_self, length(active)))
    }
    mn <- bind_rows(mn, tibble(from = self_nodes, to = self_nodes,
                               effect = as.integer(self_loop_effect),
                               strength = self_loop_strength, hill = self_loop_hill)) |>
      distinct(from, to, .keep_all = TRUE)
  }

  list(module_network = mn, self_loop_nodes = self_nodes)
}

#' Read an external topology onto a network's node ids
#'
#' Genes are mapped onto `node_ids` by a seeded random permutation. Known signs are
#' kept, unknown signs (0) are drawn 75/25 activation/repression, strengths are 1 and
#' Hill coefficients are drawn from 2 and 4 (4 for self-loops). Repressors on
#' non-self edges are capped.
#'
#' @param path tab-separated file with columns `regulator`, `target`, `sign`
#'   (1, -1 or 0 = unknown); it must contain exactly `length(node_ids)` genes.
#' @param node_ids node names of the network.
#' @param max_repressors_per_gene repressor cap, see [cap_repressors()].
#' @param seed seed.
#' @return list with `module_network`, `self_loop_nodes`, the gene-to-node map
#'   `gene_map` and `n_unknown_sign`.
#' @export
read_external_topology <- function(path, node_ids, max_repressors_per_gene = 3L, seed) {
  set.seed(seed)
  e <- utils::read.delim(path, stringsAsFactors = FALSE)
  genes <- sort(unique(c(e$regulator, e$target)))
  if (length(genes) != length(node_ids))
    stop("topology has ", length(genes), " genes but the network has ", length(node_ids), " nodes")
  map <- stats::setNames(sample(node_ids), genes)
  mn <- tibble(from = unname(map[e$regulator]), to = unname(map[e$target]), sign = e$sign)
  unknown <- mn$sign == 0
  mn$effect <- as.integer(mn$sign)
  mn$effect[unknown] <- sample(c(1L, -1L), sum(unknown), replace = TRUE, prob = c(0.75, 0.25))
  mn$strength <- 1
  mn$hill <- sample(c(2, 4), nrow(mn), replace = TRUE)
  mn$hill[mn$from == mn$to] <- 4
  mn <- mn |> select(-sign) |> distinct(from, to, .keep_all = TRUE)
  mn <- bind_rows(cap_repressors(mn[mn$from != mn$to, ], max_repressors_per_gene),
                  mn[mn$from == mn$to, ])
  list(module_network = mn, self_loop_nodes = mn$from[mn$from == mn$to],
       gene_map = tibble(gene = names(map), node = unname(map)), n_unknown_sign = sum(unknown))
}

#' Join weak fragments of a node subset into one component
#'
#' Adds activating edges from the largest fragment into the others until `members`
#' is weakly connected. Edges are only added within `members`.
#'
#' @param module_network edge table.
#' @param members nodes that must form one component.
#' @param node_ids all node names.
#' @param seed seed, or NULL.
#' @return the edge table.
#' @export
connect_fragments <- function(module_network, members, node_ids, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mn <- module_network
  if (length(members) < 2L) return(mn)
  repeat {
    g   <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
    cm  <- igraph::components(igraph::induced_subgraph(g, members), mode = "weak")
    if (cm$no <= 1L) break
    big   <- names(cm$membership)[cm$membership == which.max(cm$csize)]
    small <- names(cm$membership)[cm$membership != which.max(cm$csize)]
    tgt   <- names(cm$membership)[cm$membership == cm$membership[[small[1]]]]
    mn <- bind_rows(mn, tibble(
      from = if (length(big) > 1L) sample(big, 1) else big[1],
      to   = if (length(tgt) > 1L) sample(tgt, 1) else tgt[1],
      effect = 1L, strength = 1, hill = 2)) |>
      distinct(from, to, .keep_all = TRUE)
  }
  mn
}

#' Build the burn-in scaffolding of a network
#'
#' Per connected component: a root (the node reaching most of its component, with its
#' in-edges removed and basal expression `basal_root`), a forced activating spine of
#' `spine_len` nodes, and frontier nodes one hop downstream of the spine up to
#' `burn_frac` of the component. Self-loops lost from roots are placed on the next
#' highest out-degree nodes, and edges are added into any node the root cannot reach.
#' The burn-in never covers the whole network.
#'
#' @param module_network edge table.
#' @param node_ids node names.
#' @param spine_len nodes in each spine.
#' @param burn_frac burn-in size as a fraction of each component.
#' @param basal_root,basal_other basal expression of roots and other nodes.
#' @param component_of optional named vector assigning nodes to intended components.
#' @param n_self_loops number of self-loops to keep.
#' @param self_loop_effect,self_loop_hill,self_loop_strength self-loop attributes.
#' @param seed seed, or NULL.
#' @return list with `module_network`, dyngen `module_info`, `roots`, `spines`,
#'   `burn_ids`, `n_components`, `n_burn`, `burn_frac` and `coverage`.
#' @export
build_burn_structure <- function(module_network, node_ids, spine_len = 5L, burn_frac = 0.10,
                                 basal_root = 1.0, basal_other = 0.2, component_of = NULL,
                                 n_self_loops = NULL, self_loop_effect = 1L, self_loop_hill = 4,
                                 self_loop_strength = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mn <- module_network

  grp_list <- if (is.null(component_of)) list(node_ids) else split(node_ids, component_of[node_ids])
  for (gi in seq_along(grp_list))
    mn <- connect_fragments(mn, grp_list[[gi]], node_ids, seed = if (is.null(seed)) NULL else seed + gi)

  g <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
  comps <- if (is.null(component_of)) {
    igraph::components(g, mode = "weak")
  } else {
    mem <- component_of[node_ids]; names(mem) <- node_ids
    list(membership = mem, no = length(unique(mem)), csize = as.integer(table(mem)))
  }

  roots <- character(0); spines <- list(); burn_ids <- character(0)
  for (ci in unique(comps$membership)) {
    members <- names(comps$membership)[comps$membership == ci]
    if (length(members) < 2L) next

    reach <- vapply(members, function(v) length(igraph::subcomponent(g, v, mode = "out")), integer(1))
    root  <- names(which.max(reach))
    sp_len <- min(spine_len, length(members))
    spine  <- c(root, setdiff(sample(members), root)[seq_len(sp_len - 1L)])
    spine  <- spine[!is.na(spine)]
    roots  <- c(roots, root); spines[[as.character(ci)]] <- spine

    mn <- mn |> filter(to != root) |> filter(!(to %in% spine[-1]))
    mn <- bind_rows(tibble(from = utils::head(spine, -1), to = utils::tail(spine, -1),
                           effect = 1L, strength = 1, hill = 2), mn) |>
      distinct(from, to, .keep_all = TRUE)

    n_front <- max(0L, as.integer(floor(burn_frac * length(members))) - length(spine))
    if (n_front > 0L) {
      g2 <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
      front <- unique(unlist(lapply(spine, function(s) igraph::neighbors(g2, s, mode = "out")$name)))
      front <- setdiff(front, c(spine, burn_ids))
      if (length(front) > n_front) {
        ind   <- igraph::degree(g2, v = front, mode = "in")
        front <- names(sort(ind))[seq_len(n_front)]
      }
      burn_ids <- c(burn_ids, front)
    }
    burn_ids <- unique(c(burn_ids, spine))
  }

  if (!is.null(n_self_loops) && n_self_loops > 0L) {
    need <- n_self_loops - sum(mn$from == mn$to)
    if (need > 0L) {
      od    <- table(factor(mn$from[mn$from != mn$to], levels = node_ids))
      cands <- setdiff(node_ids[order(as.integer(od), decreasing = TRUE)],
                       c(roots, mn$from[mn$from == mn$to]))
      add   <- utils::head(cands, need)
      if (length(add))
        mn <- bind_rows(mn, tibble(from = add, to = add, effect = as.integer(self_loop_effect),
                                   strength = self_loop_strength, hill = self_loop_hill)) |>
          distinct(from, to, .keep_all = TRUE)
    }
  }

  burn_ids <- unique(burn_ids)
  stopifnot(length(burn_ids) > 0, length(burn_ids) < length(node_ids))

  module_info <- tibble(module_id = node_ids,
                        basal = ifelse(node_ids %in% roots, basal_root, basal_other),
                        burn = node_ids %in% burn_ids, independence = 1,
                        color = grDevices::rainbow(length(node_ids)))

  g3 <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
  comp_mem <- comps$membership
  n_link   <- 0L
  for (r in roots) {
    members <- names(comp_mem)[comp_mem == comp_mem[[r]]]
    repeat {
      reached   <- intersect(names(igraph::subcomponent(g3, r, mode = "out")), members)
      unreached <- setdiff(members, reached)
      if (!length(unreached) || n_link >= 5000L) break
      reach_sizes <- vapply(unreached, function(v)
        length(igraph::subcomponent(g3, v, mode = "out")), integer(1))
      tgt <- names(which.max(reach_sizes))
      src <- if (length(reached) > 1L) sample(reached, 1) else reached[1]
      mn  <- bind_rows(mn, tibble(from = src, to = tgt, effect = 1L, strength = 1, hill = 2)) |>
        distinct(from, to, .keep_all = TRUE)
      n_link <- n_link + 1L
      g3 <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
    }
  }
  reached <- unique(unlist(lapply(roots, function(r) names(igraph::subcomponent(g3, r, mode = "out")))))

  list(module_network = mn, module_info = module_info, roots = roots, spines = spines,
       burn_ids = burn_ids, n_components = comps$no, n_burn = length(burn_ids),
       burn_frac = length(burn_ids) / length(node_ids),
       coverage = length(reached) / length(node_ids))
}

#' Reconcile the edge count to an exact budget
#'
#' Removes random edges that are not spine edges or self-loops and whose removal
#' leaves both endpoints connected, or adds random edges from current regulators to
#' non-root nodes below the in-degree cap, until the network has `n_edges` edges.
#'
#' @param module_network edge table.
#' @param node_ids node names.
#' @param n_edges target edge count.
#' @param roots burn roots (never gain in-edges).
#' @param spine_edges table of spine edges (`from`, `to`), never removed.
#' @param max_in in-degree cap; NULL = none.
#' @param seed seed, or NULL.
#' @return the edge table.
#' @export
enforce_edge_budget <- function(module_network, node_ids, n_edges, roots = character(0),
                                spine_edges = NULL, max_in = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mn     <- module_network
  cap_in <- if (is.null(max_in)) length(node_ids) else max_in
  spine_key <- if (!is.null(spine_edges) && nrow(spine_edges) > 0)
    paste(spine_edges$from, spine_edges$to) else character(0)

  while (nrow(mn) > n_edges) {
    key <- paste(mn$from, mn$to)
    nonself <- mn[mn$from != mn$to, ]
    od <- table(factor(nonself$from, levels = node_ids))
    id <- table(factor(nonself$to, levels = node_ids))
    ok  <- !(key %in% spine_key) & (mn$from != mn$to) &
      (as.integer(od[mn$from]) > 1L) & (as.integer(id[mn$to]) > 1L)
    idx <- which(ok)
    if (!length(idx)) break
    drop <- sample(idx, min(length(idx), nrow(mn) - n_edges))
    mn   <- mn[-drop, , drop = FALSE]
  }

  guard <- 0L
  while (nrow(mn) < n_edges && guard < 500L) {
    guard <- guard + 1L
    need  <- n_edges - nrow(mn)
    id    <- table(factor(mn$to, levels = node_ids))
    valid_to <- setdiff(node_ids[as.integer(id) < cap_in], roots)
    od_now   <- table(factor(mn$from, levels = node_ids))
    valid_fr <- node_ids[as.integer(od_now) > 0L]
    if (!length(valid_fr)) valid_fr <- node_ids
    if (!length(valid_to) || !length(valid_fr)) break
    cand <- tibble(from = sample(valid_fr, need * 3L, replace = TRUE),
                   to   = sample(valid_to, need * 3L, replace = TRUE)) |>
      filter(from != to) |> distinct(from, to)
    cand <- anti_join(cand, mn[, c("from", "to")], by = c("from", "to"))
    if (!nrow(cand)) next
    cand <- utils::head(cand, need)
    cand$effect   <- sample(unique(mn$effect), nrow(cand), replace = TRUE)
    cand$strength <- sample(unique(mn$strength), nrow(cand), replace = TRUE)
    cand$hill     <- sample(unique(mn$hill), nrow(cand), replace = TRUE)
    mn <- bind_rows(mn, cand) |> distinct(from, to, .keep_all = TRUE)
  }
  mn
}

#' Make every node reachable from its burn root
#'
#' Adds an activating edge from a random reached node into the first unreached node
#' of the root's component until the whole component is reachable.
#'
#' @param module_network edge table.
#' @param node_ids node names.
#' @param roots burn roots.
#' @param component_of optional named vector of intended components.
#' @param seed seed, or NULL.
#' @return list with `module_network`, `n_added` and `n_unreachable_fixed`.
#' @export
repair_connectivity <- function(module_network, node_ids, roots = character(0),
                                component_of = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mn <- module_network
  before <- nrow(mn)
  n_unreach <- 0L
  if (length(roots)) {
    g <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
    comp_mem <- if (is.null(component_of)) igraph::components(g, mode = "weak")$membership else {
      m <- component_of[node_ids]; names(m) <- node_ids; m
    }
    guard <- 0L
    for (r in roots) {
      members <- names(comp_mem)[comp_mem == comp_mem[[r]]]
      repeat {
        reached   <- intersect(names(igraph::subcomponent(g, r, mode = "out")), members)
        unreached <- setdiff(members, reached)
        if (!length(unreached) || guard >= 500L) break
        guard <- guard + 1L; n_unreach <- n_unreach + 1L
        src <- if (length(reached) > 1L) sample(reached, 1) else reached[1]
        mn  <- bind_rows(mn, tibble(from = src, to = unreached[1], effect = 1L, strength = 1,
                                    hill = 2)) |>
          distinct(from, to, .keep_all = TRUE)
        g <- graph_from_data_frame(mn[, c("from", "to")], vertices = node_ids, directed = TRUE)
      }
    }
  }
  list(module_network = mn, n_added = nrow(mn) - before, n_unreachable_fixed = n_unreach)
}

#' Summary statistics of a network
#'
#' @param module_network edge table.
#' @param node_ids node names.
#' @param burn_ids burn-in nodes.
#' @return one-row data frame: edges, sparsity with and without the diagonal,
#'   self-loops, out-degree mean / median / max / CV / Gini, leaf and source shares,
#'   isolated nodes, maximum in-degree, weak components and burn-in size.
#' @export
topology_stats <- function(module_network, node_ids, burn_ids = character(0)) {
  od <- as.numeric(table(factor(module_network$from, levels = node_ids)))
  id <- as.numeric(table(factor(module_network$to, levels = node_ids)))
  g  <- graph_from_data_frame(module_network[, c("from", "to")], vertices = node_ids, directed = TRUE)
  gini <- function(x) {
    x <- sort(x); k <- length(x)
    if (sum(x) == 0) return(NA_real_)
    sum((2 * seq_len(k) - k - 1) * x) / (k * sum(x))
  }
  n  <- length(node_ids)
  ne <- nrow(module_network)
  data.frame(n_nodes = n, n_edges = ne, sparsity = 1 - ne / n^2,
             sparsity_nodiag = 1 - ne / (n * (n - 1)),
             n_self_loops = sum(module_network$from == module_network$to),
             max_out = max(od), mean_out = mean(od), cv_out = stats::sd(od) / mean(od),
             gini_out = gini(od), median_out = stats::median(od),
             frac_zero_out = mean(od == 0), frac_zero_in = mean(id == 0),
             n_isolated = sum(od == 0 & id == 0), max_in = max(id),
             n_components = igraph::components(g, mode = "weak")$no,
             n_burn = length(burn_ids), burn_frac = length(burn_ids) / n)
}
