test_that("every bundled network loads with matching designs", {
  for (nm in bundled_networks()) {
    net <- load_network(nm)
    tf <- model_tf_ids(net$model)
    expect_equal(length(tf), 222L)
    expect_equal(colnames(net$design_random), tf)
    expect_equal(colnames(net$design_hub), tf)
    expect_equal(dim(net$grn_truth), c(222L, 222L))
    expect_gte(length(net$seeds), 300L)
  }
})

test_that("bundled dataset configs give the published screen sizes", {
  cfgs <- bundled_configs("datasets")
  expect_length(cfgs, 22L)
  for (nm in names(cfgs)) {
    cfg <- utils::modifyList(default_dataset_config(), read_config(nm))
    d <- screen_design(load_network(cfg$network), cfg)
    n <- if (nm == "legacy_knockdown") 76L else if (cfg$perturbation == "knockdown") 380L else 228L
    expect_equal(nrow(d$onehot), n, info = nm)
  }
})

test_that("network configs name a bundled network or topology", {
  for (nm in names(bundled_configs("networks"))) {
    cfg <- read_config(nm)
    expect_identical(cfg$type, "network")
    expect_true(nm %in% bundled_networks())
    if (cfg$topology == "external") expect_true(file.exists(topology_path(cfg$topology_file)))
  }
})
