context("BPOT — evd_wrappers helpers")

test_that("normalize_eta works", {
  expect_equal(normalize_eta(0.5), c(0.5, 0.5))
  expect_equal(normalize_eta(c(0.3, 0.7)), c(0.3, 0.7))
})

test_that("tail_adjust returns values in [0,1] for all tail types", {
  x1 <- 1; x2 <- 2; v <- x1 + x2
  p2 <- tail_adjust(v, x1, x2, 2)
  p3 <- tail_adjust(v, x1, x2, 3)
  p4 <- tail_adjust(v, x1, x2, 4)
  expect_true(p2 >= 0 && p2 <= 1)
  expect_true(p3 >= 0 && p3 <= 1)
  expect_true(p4 >= 0 && p4 <= 1)
})

test_that("tail_adjust type-2 is positive", {
  x1 <- 0.8; x2 <- 1.2; v <- x1 + x2
  p2 <- tail_adjust(v, x1, x2, 2)
  expect_gt(p2, 0)
  expect_lt(p2, 1)
})

test_that("mtransform_gp_mk2 transforms to exp margin", {
  mar <- c(1.0, -0.1)
  thres <- -2
  eta <- 0.15
  q <- seq(-1.5, -0.5, length.out = 5)
  out <- mtransform_gp_mk2(q, p = mar, thres = thres, eta = eta, margin = "exp")
  expect_type(out, "double")
  expect_true(all(out >= 0))
})

test_that("mtransform_gp_mk2 transforms to uniform margin in [0,1]", {
  mar <- c(1.0, -0.1)
  thres <- -2
  eta <- 0.15
  q <- seq(-1.5, -0.5, length.out = 5)
  out <- mtransform_gp_mk2(q, p = mar, thres = thres, eta = eta, margin = "uniform")
  expect_true(all(out >= 0 & out <= 1))
})

test_that("mtransform_gp_mk2 transforms to frechet margin", {
  mar <- c(1.0, -0.1)
  thres <- -2
  eta <- 0.15
  q <- seq(-1.5, -0.5, length.out = 5)
  out <- mtransform_gp_mk2(q, p = mar, thres = thres, eta = eta, margin = "frechet")
  expect_type(out, "double")
  expect_true(all(out > 0))
})

test_that("mtransform_gp_mk2 errors on values below threshold", {
  mar <- c(1.0, -0.1)
  thres <- -1
  eta <- 0.15
  expect_error(mtransform_gp_mk2(-2, p = mar, thres = thres, eta = eta),
               "input below thresholds")
})

test_that("mtransform_gp_mk2 errors on invalid margin", {
  mar <- c(1.0, -0.1)
  thres <- -3
  eta <- 0.15
  expect_error(mtransform_gp_mk2(-1.5, p = mar, thres = thres, eta = eta,
                                 margin = "invalid"),
               "invalid margin type")
})


context("BPOT — bivariate CDF (pb_tv*) wrappers")

test_that("pb_tvevd dispatches to logistic model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "log", dep = 0.6,
                mar1 = mar1_synth, mar2 = mar2_synth, tail_type = 2,
                thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd dispatches to HR model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "hr", dep = 1.5,
                mar1 = mar1_synth, mar2 = mar2_synth, tail_type = 2,
                thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd dispatches to negative logistic model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "neglog",
                dep = 1.5, mar1 = mar1_synth, mar2 = mar2_synth,
                tail_type = 2, thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd dispatches to asymmetric logistic model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "alog", dep = 0.6,
                asy = c(0.8, 0.9), mar1 = mar1_synth, mar2 = mar2_synth,
                tail_type = 2, thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd dispatches to bilogistic model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "bilog",
                alpha = 0.7, beta = 0.8, mar1 = mar1_synth,
                mar2 = mar2_synth, tail_type = 2, thres = thres,
                eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd dispatches to negative bilogistic model", {
  thres <- c(u_synth, v_synth)
  p <- pb_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "negbilog",
                alpha = 1.2, beta = 0.8, mar1 = mar1_synth,
                mar2 = mar2_synth, tail_type = 2, thres = thres,
                eta = c(eta_synth, eta_synth))
  expect_true(is.finite(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("pb_tvevd rejects invalid models", {
  expect_error(pb_tvevd(q1 = 0, q2 = 0, model = "xyz"),
               "should be one of")
})


context("BPOT — bivariate density (db_tv*) wrappers")

test_that("db_tvevd returns finite value for logistic model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "log",
                dep = 0.6, mar1 = mar1_synth, mar2 = mar2_synth,
                thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(d))
})

test_that("db_tvevd returns finite value for HR model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "hr",
                dep = 1.5, mar1 = mar1_synth, mar2 = mar2_synth,
                thres = thres, eta = c(eta_synth, eta_synth),
                margin = "exp")
  expect_true(is.finite(d))
})

test_that("db_tvevd returns finite value for negative logistic model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "neglog",
                dep = 1.5, mar1 = mar1_synth, mar2 = mar2_synth,
                thres = thres, eta = c(eta_synth, eta_synth))
  expect_true(is.finite(d))
})

test_that("db_tvevd returns finite value for asymmetric logistic model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "alog",
                dep = 0.6, asy = c(0.8, 0.9), mar1 = mar1_synth,
                mar2 = mar2_synth, thres = thres,
                eta = c(eta_synth, eta_synth))
  expect_true(is.finite(d))
})

test_that("db_tvevd returns finite value for bilogistic model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "bilog",
                alpha = 0.7, beta = 0.8, mar1 = mar1_synth,
                mar2 = mar2_synth, thres = thres,
                eta = c(eta_synth, eta_synth))
  expect_true(is.finite(d))
})

test_that("db_tvevd returns finite value for negative bilogistic model", {
  thres <- c(u_synth, v_synth)
  d <- db_tvevd(q1 = u_synth + 1, q2 = v_synth + 5, model = "negbilog",
                alpha = 1.2, beta = 0.8, mar1 = mar1_synth,
                mar2 = mar2_synth, thres = thres,
                eta = c(eta_synth, eta_synth))
  expect_true(is.finite(d))
})


context("BPOT — a_bp_approx")

test_that("a_bp_approx returns finite values for ord='0'", {
  beta <- c(1, 1.1, 0.9, 0.8, 0.85)
  vals <- sapply(seq(0.1, 0.9, length.out = 5), a_bp_approx, A_bp = beta,
                 ord = "0")
  expect_true(all(is.finite(as.numeric(vals))))
})


context("BPOT — integrate_safe")

test_that("integrate_safe returns finite value for simple integrand", {
  f <- function(x) dnorm(x)
  val <- integrate_safe(f, lower = -1, upper = 1)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})

test_that("integrate_safe returns NA on non-finite integrand", {
  f <- function(x) ifelse(x > 0.5, Inf, 1)
  val <- integrate_safe(f, lower = 0, upper = 1)
  expect_true(is.na(val))
})


context("BPOT — bpot_dep_change")

test_that("bpot_dep_change modifies logistic dependence parameter", {
  M <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  M_mod <- bpot_dep_change(0.5, M$M)
  expect_equivalent(M_mod$estimate[5], 0.5)
  expect_equivalent(M_mod$estimate[1:4], M$M$estimate[1:4])
})

test_that("bpot_dep_change modifies HR dependence parameter", {
  M <- run_single_bpot(synth_bpot, model = "hr",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  M_mod <- bpot_dep_change(1.2, M$M)
  expect_equivalent(M_mod$estimate[5], 1.2)
})


context("BPOT — run_single_bpot + create_result_bpot")

test_that("run_single_bpot returns correct list structure", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  expect_type(r, "list")
  expect_true(all(c("M", "pcrash", "plot_df", "ss", "xcrash") %in% names(r)))
  expect_s3_class(r$M, "bvpot")
  expect_type(r$pcrash, "double")
  expect_s3_class(r$plot_df, "data.frame")
})

test_that("run_single_bpot works for log model", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  expect_equal(r$M$model, "log")
  expect_true(r$pcrash > 0 && r$pcrash < 1)
})

test_that("run_single_bpot works for HR model", {
  r <- run_single_bpot(synth_bpot, model = "hr",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  expect_equal(r$M$model, "hr")
  expect_true(r$pcrash > 0 && r$pcrash < 1)
})

test_that("create_result_bpot combines single-model results into a list", {
  r1 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  r2 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  res <- create_result_bpot(CN = r1, SE = r2)
  expect_type(res, "list")
  expect_length(res, 2)
  expect_equal(names(res), c("CN", "SE"))
  expect_true(all(c("M", "pcrash", "plot_df") %in% names(res[[1]])))
})


context("BPOT — create_plot_df")

test_that("create_plot_df returns valid data frame for log model", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  df <- create_plot_df(r$ss, x = r$xcrash, model = r$M, px = r$pcrash)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("speed", "JointP", "ConditionP", "ConditionalD") %in%
                  names(df)))
  expect_true(all(df$ConditionalD >= 0, na.rm = TRUE))
})


context("BPOT — call_c_bivariate and normalize_c_bivariate")

test_that("call_c_bivariate returns finite positive value", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  y_val <- max(v_synth + 5, 20)
  val <- call_c_bivariate(y = y_val, x = x0_synth, ev_model = r$M)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})

test_that("normalize_c_bivariate returns finite positive value", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  val <- normalize_c_bivariate(x = x0_synth, ev_model = r$M)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})


context("BPOT — injury probability")

test_that("injury_from_c_bivariate returns finite positive value", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  ip <- injury_from_c_bivariate(r$plot_df$speed, ev_model = r$M,
                                severity = pis0_test, x0 = x0_synth,
                                px = r$pcrash)
  expect_true(is.finite(ip))
  expect_gt(ip, 0)
})

test_that("injury_from_c_bivariate_e returns finite positive value", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  ip <- injury_from_c_bivariate_e(r$plot_df$speed, ev_model = r$M,
                                  severity = pis1_test, x0 = x0_synth,
                                  px = r$pcrash, N = 50, age_mean = 40)
  expect_true(is.finite(ip))
  expect_gt(ip, 0)
})


context("BPOT — summarise_bpot")

test_that("summarise_bpot runs without error", {
  r1 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  r2 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  res <- create_result_bpot(r1, r2)
  expect_output(
    summarise_bpot(res, model_names = c("log", "log")),
    "BPOT Model Summary"
  )
})


context("BPOT — theoretical density")

test_that("bpot_theoretical_density returns data frame with origin column", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  tf <- td_tempfile("bpot")
  td <- bpot_theoretical_density(c(0.6, 0.8), r, filename = basename(tf))
  expect_s3_class(td, "data.frame")
  expect_true("origin" %in% names(td))
})


context("BPOT — plot functions")

test_that("plot_crash_severity returns density and distribution plots", {
  r1 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  r2 <- run_single_bpot(synth_bpot, model = "log",
                        thres = c(u_synth, v_synth), xcrash = x0_synth)
  p <- plot_crash_severity(list(r1$plot_df, r2$plot_df), injury_df_test)
  expect_type(p, "list")
  expect_s3_class(p$density, "ggplot")
  expect_s3_class(p$distribution, "ggplot")

  p_d <- plot_crash_severity(list(r1$plot_df, r2$plot_df), injury_df_test,
                             type = "density")
  expect_s3_class(p_d, "ggplot")

  p_F <- plot_crash_severity(list(r1$plot_df, r2$plot_df), injury_df_test,
                             type = "distribution")
  expect_s3_class(p_F, "ggplot")
})

test_that("plot_crash_severity_single returns a ggplot object", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  p <- plot_crash_severity_single(r$plot_df, injury_df_test)
  expect_s3_class(p, "ggplot")
})

test_that("plot_theoretical_density returns a ggplot object", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  tf <- td_tempfile("bpot_plot")
  td <- bpot_theoretical_density(c(0.6, 0.8), r, filename = basename(tf))
  p <- plot_theoretical_density(td)
  expect_s3_class(p, "ggplot")
})


context("BPOT — gof_bpot")

test_that("gof_bpot returns htest object", {
  r <- run_single_bpot(synth_bpot, model = "log",
                       thres = c(u_synth, v_synth), xcrash = x0_synth)
  gf <- gof_bpot(r$M, synth_bpot)
  expect_s3_class(gf, "htest")
  expect_true("p.value" %in% names(gf))
})


context("BPOT — pickands_nonpar")

test_that("pickands_nonpar returns list with beta coefficients", {
  eta <- nrow(subset(synth_bpot, prox > u_synth &
                     Speed > v_synth)) / nrow(synth_bpot)
  thres <- c(u_synth, v_synth)
  res <- pickands_nonpar(dat = synth_bpot, mar1 = mar1_synth,
                         mar2 = mar2_synth,
                         thres = thres, eta = eta,
                         est = "cfg", CI = FALSE, d = 2, k = 10,
                         N = 100, ifplot = FALSE)
  expect_type(res, "list")
  expect_true("beta" %in% names(res))
})
