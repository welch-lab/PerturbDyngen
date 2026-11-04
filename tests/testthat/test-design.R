test_that("hub design perturbs the activating regulators of one hub per experiment", {
  g <- paste0("G", 1:10)
  adj <- matrix(0, 10, 10, dimnames = list(g, g))
  adj[1:5, "G10"] <- 1
  adj[6:7, "G9"] <- 1
  adj[8, "G9"] <- -1
  h <- hub_design(adj, max_nodes = 3L, top_k = 2L, min_upstream = 2L)
  expect_equal(attr(h, "hub_of_experiment"), c("G10", "G10", "G9"))
  expect_equal(unname(rowSums(h)), c(3, 2, 2))
  expect_equal(h[3, "G8"], 0L)
})

test_that("fair_ordered covers every gene", {
  fo <- fair_ordered(20, 8, max_nodes_per_experiment = 10, seed = 1)
  expect_equal(dim(fo$experiment_matrix), c(8L, 20L))
  expect_true(all(colSums(fo$experiment_matrix) == fo$l))
})
