context("CI — exposure argument")

test_that("all CI methods expose an exposure argument", {
  fns <- c("ci_delta", "ci_prolik", "ci_bayes", "ci_sim",
           "ci_bpot_sim_joint", "ci_bpot_sim_marginal", "ci_cbpot_sim")
  for (nm in fns) {
    expect_true("exposure" %in% names(formals(get(nm))), info = nm)
  }
})
