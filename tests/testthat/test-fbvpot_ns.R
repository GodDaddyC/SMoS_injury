context("BPOT — fbvpot_ns (non-stationary scale)")

make_test_df <- function(n = 150, seed = 1) {
  set.seed(seed)
  x <- cbind(rexp(n, .3), rgamma(n, 3, .08))
  data.frame(
    x1 = x[, 1],
    x2 = x[, 2],
    v_conflict = seq(-1, 1, length.out = n),
    movement = factor(sample(0:3, n, TRUE), levels = 0:3, ordered = TRUE)
  )
}

test_that("fbvpot_ns matches evd::fbvpot in the stationary case", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  x <- as.matrix(df[, c("x1", "x2")])

  fb <- fbvpot(x, threshold = u, model = "log")
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u)

  expect_equal(as.numeric(ns$estimate), as.numeric(fb$estimate), tolerance = 1e-4)
  expect_equal(ns$deviance, fb$deviance, tolerance = 1e-8)
  expect_equal(ns$nat, fb$nat)
  expect_equal(ns$threshold, fb$threshold)
})

test_that("fbvpot_ns returns evd-style list with design fields", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict", nsscale2 = "v_conflict")

  expect_type(ns, "list")
  expect_true(all(c("estimate", "std.err", "deviance", "var.cov",
                    "threshold", "nat", "n", "model", "nsscale1",
                    "nsscale2", "design1", "design2", "scale1_terms",
                    "scale2_terms", "resp") %in% names(ns)))
  expect_equal(names(ns$estimate),
               c("scale1_0", "scale1_1", "shape1",
                 "scale2_0", "scale2_1", "shape2", "dep"))
  expect_equal(ns$n, 150)
  expect_equal(ns$model, "log")
  expect_equal(ns$resp, c("x1", "x2"))
})

test_that("fbvpot_ns handles multiple numeric covariates", {
  df <- make_test_df()
  df$v2 <- df$v_conflict^2
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = c("v_conflict", "v2"), std.err = FALSE)
  expect_equal(names(ns$estimate),
               c("scale1_0", "scale1_1", "scale1_2", "shape1",
                 "scale2_0", "shape2", "dep"))
  expect_equal(ns$scale1_terms$covariate, c("(Intercept)", "v_conflict", "v2"))
})

test_that("fbvpot_ns handles factor covariates with treatment contrasts", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "movement", std.err = FALSE)
  expect_equal(names(ns$estimate),
               c("scale1_0", "scale1_1", "scale1_2", "scale1_3", "shape1",
                 "scale2_0", "shape2", "dep"))
  expect_equal(ns$scale1_terms$covariate,
               c("(Intercept)", "movement", "movement", "movement"))
  expect_equal(ns$scale1_terms$level, c(NA, "1", "2", "3"))

  # treatment-coded design matches manual indicators (reference level 0 dropped)
  mm <- cbind(1,
              as.numeric(df$movement == "1"),
              as.numeric(df$movement == "2"),
              as.numeric(df$movement == "3"))
  expect_equal(unname(ns$design1), unname(mm))
})

test_that("fbvpot_ns handles mixed numeric + factor covariates per margin", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = c("v_conflict", "movement"),
                  nsscale2 = "movement", std.err = FALSE)
  expect_equal(names(ns$estimate),
               c("scale1_0", "scale1_1", "scale1_2", "scale1_3", "scale1_4",
                 "shape1",
                 "scale2_0", "scale2_1", "scale2_2", "scale2_3",
                 "shape2", "dep"))
  expect_equal(ns$scale1_terms$covariate,
               c("(Intercept)", "v_conflict", "movement", "movement", "movement"))
  expect_true(is.na(ns$scale1_terms$level[2]))
})

test_that("fbvpot_ns handles per-margin covariates independently", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns1 <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                   nsscale1 = "v_conflict")
  expect_equal(names(ns1$estimate),
               c("scale1_0", "scale1_1", "shape1",
                 "scale2_0", "shape2", "dep"))
  expect_null(ns1$nsscale2)
})

test_that("fbvpot_ns keeps scale positive under identity link", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  sc <- drop(ns$design1 %*% ns$estimate[ns$scale1_terms$param])
  expect_true(all(sc > 0))
})

test_that("fbvpot_ns validates inputs", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))

  expect_error(fbvpot_ns(as.matrix(df), resp = c("x1", "x2"), threshold = u),
               "data.frame")
  expect_error(fbvpot_ns(df, resp = "x1", threshold = u),
               "length two")
  expect_error(fbvpot_ns(df, resp = c("x1", "nope"), threshold = u),
               "not found in data")
  expect_error(fbvpot_ns(df, resp = c("x1", "x2")),
               "threshold")
  expect_error(fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                         nsscale1 = "nope"),
               "not found in data")
})

test_that("fbvpot_ns rejects non-logistic models", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  expect_error(fbvpot_ns(df, resp = c("x1", "x2"), threshold = u, model = "hr"),
               "should be one of|\"log\"")
})

test_that("fbvpot_ns recovers non-stationary scale signal", {
  set.seed(42)
  n <- 300
  cov <- seq(-1, 1, length.out = n)
  u1 <- 3
  u2 <- 50
  s1 <- exp(0.5 + 1.0 * cov)
  shape1 <- 0.1
  x1 <- numeric(n)
  for (i in seq_len(n)) {
    if (runif(1) < 0.8) {
      x1[i] <- runif(1, 0, u1)
    } else {
      x1[i] <- u1 + rexp(1, rate = (1 - shape1) / s1[i])
    }
  }
  x2 <- ifelse(runif(n) < 0.8, runif(n, 0, u2), u2 + rexp(n, 0.05))
  df <- data.frame(x1 = x1, x2 = x2, v_conflict = cov)

  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = c(u1, u2),
                  nsscale1 = "v_conflict")
  expect_gt(ns$estimate["scale1_1"], 0.5)
})


context("BPOT — fbvpot_ns methods (fitted, confint, plot, print)")

test_that("fbvpot_ns returns a classed object", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  expect_s3_class(ns, "fbvpot_ns")
})

test_that("fitted.fbvpot_ns returns per-observation scale per margin", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  f <- fitted(ns)
  expect_type(f, "list")
  expect_length(f, 2)
  expect_equal(length(f[[1]]), nrow(df))
  expect_equal(length(f[[2]]), nrow(df))
  expect_true(all(f[[1]] > 0))
})

test_that("confint.fbvpot_ns returns normal-approximation intervals", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  ci <- confint(ns)
  expect_true(all(dim(ci) == c(length(ns$estimate), 2)))
  expect_true(all(ci[, 1] <= ci[, 2]))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))

  ci2 <- confint(ns, parm = c("scale1_1", "dep"))
  expect_equal(rownames(ci2), c("scale1_1", "dep"))

  # centre is the estimate
  mid <- (ci2[, 1] + ci2[, 2]) / 2
  expect_equal(mid, ns$estimate[c("scale1_1", "dep")], tolerance = 1e-6)
})

test_that("confint.fbvpot_ns fails when std.err is unavailable", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict", std.err = FALSE)
  expect_error(confint(ns), "std.err")
})

test_that("confint.fbvpot_ns rejects unknown parameters", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u)
  expect_error(confint(ns, parm = "nope"), "unknown parameters")
})

test_that("plot.fbvpot_ns returns a ggplot for each num", {
  df <- make_test_df(n = 400)
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = c("v_conflict", "movement"),
                  nsscale2 = "v_conflict")

  p1 <- plot(ns, num = 1)
  expect_s3_class(p1, "ggplot")
  p2 <- plot(ns, num = 2)
  expect_s3_class(p2, "ggplot")
  p3 <- plot(ns, num = 3)
  expect_s3_class(p3, "ggplot")
  p4 <- plot(ns, num = 4)
  expect_s3_class(p4, "ggplot")
})

test_that("plot.fbvpot_ns errors on invalid num", {
  df <- make_test_df(n = 400)
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  expect_error(plot(ns, num = 5), "subset of 1:4")
})

test_that("print.fbvpot_ns runs without error", {
  df <- make_test_df()
  u <- c(quantile(df$x1, .8), quantile(df$x2, .8))
  ns <- fbvpot_ns(df, resp = c("x1", "x2"), threshold = u,
                  nsscale1 = "v_conflict")
  expect_output(print(ns), "Bivariate POT")
})
