context("CBPOT — q_distr_param")

test_that("q_distr_param returns mvdc for type 1", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  mv <- q_distr_param(cop, mar1 = pot_synth, mar2 = conseq_synth, type = 1)
  expect_s4_class(mv, "mvdc")
})

test_that("q_distr_param returns mvdc for type 2", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  mv <- q_distr_param(cop, mar1 = pot_synth, mar2 = conseq_synth, type = 2)
  expect_s4_class(mv, "mvdc")
})

test_that("q_distr_param returns mvdc for type 4", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  mv <- q_distr_param(cop, mar1 = pot_synth, mar2 = conseq_synth, type = 4)
  expect_s4_class(mv, "mvdc")
})


context("CBPOT — c_q_bivariate and normalize_c_q_bivariate")

test_that("c_q_bivariate returns finite positive value", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  val <- c_q_bivariate(y = 0.3, x = 0.8, model = cop)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})

test_that("c_q_bivariate returns zero at upper bound", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  val <- c_q_bivariate(y = 0.5, x = 1, model = cop)
  expect_equal(val, 0, tolerance = 1e-6)
})

test_that("normalize_c_q_bivariate returns finite positive value", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  val <- normalize_c_q_bivariate(x = 0.8, model = cop, lb = 0)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})


context("CBPOT — run_single_cbpot + create_result_cbpot")

test_that("run_single_cbpot returns correct list structure (gumbel)", {
  r <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  expect_type(r, "list")
  expect_true(all(c("Cop", "CM", "CM0", "x0_un", "pcrash", "pu",
                    "p2", "s_un", "plot_df_q") %in% names(r)))
  expect_s4_class(r$Cop, "fitCopula")
  expect_s4_class(r$CM, "mvdc")
  expect_type(r$pcrash, "double")
  expect_gt(r$pcrash, 0)
  expect_lt(r$pcrash, 1)
})

test_that("run_single_cbpot works for clayton copula", {
  r <- run_single_cbpot(cop_dat_synth, copula = claytonCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth)
  expect_s4_class(r$Cop, "fitCopula")
  expect_true(r$pcrash > 0 && r$pcrash < 1)
})

test_that("run_single_cbpot works for normal copula", {
  r <- run_single_cbpot(cop_dat_synth,
                        copula = normalCopula(param = -0.1, dispstr = "ex"),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  expect_s4_class(r$Cop, "fitCopula")
  expect_true(r$pcrash > 0 && r$pcrash < 1)
})

test_that("create_result_cbpot combines single-model results into a list", {
  r1 <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                         p2 = conseq_synth, qcrash = qcrash_synth,
                         pot = pot_synth, pu = pu_synth, method = "mpl")
  r2 <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                         p2 = conseq_synth, qcrash = qcrash_synth,
                         pot = pot_synth, pu = pu_synth, method = "mpl")
  res <- create_result_cbpot(CN = r1, SE = r2)
  expect_type(res, "list")
  expect_length(res, 2)
  expect_equal(names(res), c("CN", "SE"))
  expect_true(all(c("Cop", "CM", "CM0", "plot_df_q") %in% names(res[[1]])))
})


context("CBPOT — create_plot_df_q")

test_that("create_plot_df_q returns valid data frame", {
  r <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  df <- r$plot_df_q
  expect_s3_class(df, "data.frame")
  expect_true(all(c("speed_u", "speed", "JointP", "ConditionP",
                    "ConditionalD") %in% names(df)))
  expect_true(all(df$ConditionalD >= 0, na.rm = TRUE))
})

test_that("create_plot_df_q length matches s_un", {
  r <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  expect_equal(nrow(r$plot_df_q), length(s_un_synth))
})


context("CBPOT — injury probability")

test_that("injury_from_c_q_bivariate returns probability", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  ip <- injury_from_c_q_bivariate(cop, pis0_test, x0 = 0.8, px = 0.15,
                                  p2 = conseq_synth, pu = 0.15)
  expect_true(is.finite(ip))
  expect_gt(ip, 0)
})

test_that("injury_from_c_q_bivariate_e returns finite value", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  ip <- injury_from_c_q_bivariate_e(cop, pis1_test, x0 = 0.8, px = 0.15,
                                    p2 = conseq_synth, pu = 0.15,
                                    N = 100, age_mean = 40)
  expect_true(is.finite(ip))
})


context("CBPOT — summarise_cbpot")

test_that("summarise_cbpot runs without error", {
  r1 <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                         p2 = conseq_synth, qcrash = qcrash_synth,
                         pot = pot_synth, pu = pu_synth, method = "mpl")
  r2 <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                         p2 = conseq_synth, qcrash = qcrash_synth,
                         pot = pot_synth, pu = pu_synth, method = "mpl")
  res <- create_result_cbpot(r1, r2)
  expect_output(
    summarise_cbpot(res, pot = list(pot_synth, pot_synth),
                    model_names = c("gumbel", "gumbel")),
    "CBPOT Model Summary"
  )
})


context("CBPOT — theoretical density")

test_that("cbpot_theoretical_density returns data frame with origin column", {
  r <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  tf <- td_tempfile("cbpot")
  td <- cbpot_theoretical_density(c(1.3, 1.8), "gumbel", r,
                                  filename = basename(tf))
  expect_s3_class(td, "data.frame")
  expect_true("origin" %in% names(td))
})

context("CBPOT — theoretical severe probability")

test_that("cbpot_theoretical_severe returns conditional probabilities", {
  r <- run_single_cbpot(cop_dat_synth, copula = gumbelCopula(),
                        p2 = conseq_synth, qcrash = qcrash_synth,
                        pot = pot_synth, pu = pu_synth, method = "mpl")
  params <- c(1.2, 1.8)
  out <- cbpot_theoretical_severe(params, "gumbel", r,
                                  severity_boundary = 40)
  speed_tail <- pgamma(40, shape = r$p2$estimate[1],
                       rate = r$p2$estimate[2], lower.tail = FALSE)
  expected <- vapply(params, function(param) {
    copula <- gumbelCopula(param = param)
    pCopula(c(r$qcrash, speed_tail), copula)
  }, numeric(1))

  expect_type(out, "double")
  expect_length(out, length(params))
  expect_equal(out, expected)
  expect_true(all(out >= 0 & out <= 1))
})


context("CBPOT — nonparametric helpers")

test_that("c_q_bivariate_nonpar returns finite value", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  u <- rCopula(200, cop)
  kd <- kdecop(u)
  val <- c_q_bivariate_nonpar(v = 0.3, u = 0.8, model = kd)
  expect_true(is.finite(val))
})

test_that("create_plot_df_q_nonpar returns valid data frame", {
  cop <- gumbelCopula(param = 1.5, dim = 2)
  u <- rCopula(200, cop)
  kd <- kdecop(u)
  df <- create_plot_df_q_nonpar(s_un_synth, x = 0.8, model = kd,
                                px = 0.15, p2 = conseq_synth)
  expect_s3_class(df, "data.frame")
  expect_true("ConditionalD" %in% names(df))
})
