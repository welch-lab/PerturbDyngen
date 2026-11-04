test_that("sigmoid inverse round-trips", {
  kd <- c(0.005, 0.2, 0.5, 0.8)
  expect_equal(sigmoid_kd(sigmoid_inverse(kd)), kd)
})

test_that("guide UMIs are reproducible and bounded by capture", {
  a <- simulate_guide_umis(50, c("G1", "G2"), kd_level = 0.2, seed = 7)
  b <- simulate_guide_umis(50, c("G1", "G2"), kd_level = 0.2, seed = 7)
  expect_identical(a, b)
  expect_true(all(a$obs_umis <= a$true_umis))
  expect_true(all(is.na(simulate_guide_umis(5, "G1", 0.2, u_target = 200, seed = 1)$kd_continuous)))
})

test_that("RNA UMIs never exceed true counts", {
  mol <- c(paste0("mol_premrna_G", 1:3), paste0("mol_mrna_G", 1:3))
  mat <- array(rpois(6 * 2 * 4, 30), c(6, 2, 4))
  set.seed(1)
  r <- simulate_rna_umis(mat, mol, 2, 0.5, 0.1, 3, rep(1, 4))
  expect_equal(dim(r$obs_spliced), c(4L, 3L))
  expect_true(all(r$obs_spliced <= t(mat[4:6, 2, ])))
  expect_true(all(r$obs_unspliced <= t(mat[1:3, 2, ])))
})
