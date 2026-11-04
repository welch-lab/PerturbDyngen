test_that("a paired cell is identical to wild type before the perturbation", {
  net <- load_network("legacy")
  sys <- simulation_setup(net$model, pre_window = 0.05)
  genes <- colnames(net$design_random)[net$design_random[1, ] == 1]
  res <- simulate_experiment(genes, 0.005, sys, net$seeds[1])
  pre <- res$snap_times < sys$pert_time
  expect_equal(sum(pre), 4L)
  expect_identical(res$mat_wt[, pre, 1], res$mat_kd[, pre, 1])
  mrna <- paste0("mol_mrna_", genes)
  last <- length(res$snap_times)
  expect_lt(sum(res$mat_kd[match(mrna, res$mol_ids), last, 1]),
            sum(res$mat_wt[match(mrna, res$mol_ids), last, 1]))
})
