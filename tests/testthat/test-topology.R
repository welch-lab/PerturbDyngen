test_that("sampled networks hit the edge budget and are reachable from the root", {
  nodes <- paste0("G", 1:60)
  topo <- sample_topology(nodes, n_edges = 150, sdlog = 1, frac_zero_out = 0.2,
                          frac_self_loops = 0.05, seed = 3)
  burn <- build_burn_structure(topo$module_network, nodes,
                               n_self_loops = length(topo$self_loop_nodes), seed = 3)
  spine <- dplyr::bind_rows(lapply(burn$spines, function(s)
    tibble::tibble(from = utils::head(s, -1), to = utils::tail(s, -1))))
  mn <- enforce_edge_budget(burn$module_network, nodes, 150, roots = burn$roots, spine_edges = spine, seed = 3)
  mn <- repair_connectivity(mn, nodes, roots = burn$roots, seed = 3)$module_network
  st <- topology_stats(mn, nodes, burn$burn_ids)
  expect_lte(abs(st$n_edges - 150), 3)
  expect_equal(st$n_components, 1L)
  expect_equal(burn$coverage, 1)
  expect_true(all(mn$to != burn$roots | mn$from == mn$to))
})

test_that("repressor cap holds", {
  mn <- tibble::tibble(from = paste0("G", 1:6), to = "T", effect = -1L, strength = 1, hill = 2)
  expect_equal(sum(cap_repressors(mn, 3L, seed = 1)$effect == -1L), 3L)
})
