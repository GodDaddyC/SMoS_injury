# wrapper function for fevd object in extremes

pit_gp <- function(x, threshold, scale, shape,lower.tail=FALSE) {
  # vetorizatio for evaluation of non-stationary parameters
  
  n <- length(x)
  if (length(scale) == 1) scale <- rep(scale, n)
  if (length(shape) == 1) shape <- rep(shape, n)
  mapply(pevd, x, threshold = threshold, scale = scale, shape = shape,
         MoreArgs = list(type = "GP",lower.tail=lower.tail))
}

pit_gev <- function(x, loc, scale, shape,lower.tail=TRUE) {
  # vetorizatio for evaluation of non-stationary parameters
  
  n <- length(x)
  if (length(loc) == 1) loc <- rep(loc, n)
  if (length(scale) == 1) scale <- rep(scale, n)
  if (length(shape) == 1) shape <- rep(shape, n)
  mapply(pevd, x,loc = loc,scale = scale, shape = shape,
         MoreArgs = list(type = "GEV",lower.tail=lower.tail))
}

build_ns_param <- function(term_names, cov_data, par, prefix) {
  n <- nrow(cov_data)
  intercept <- par[[paste0(prefix, "0")]]

  if (is.null(term_names) || length(term_names) == 0)
    return(rep(intercept, n))

  result <- rep(intercept, n)
  pidx <- 1L

  for (i in seq_along(term_names)) {
    col_name <- term_names[i]

    if (!col_name %in% names(cov_data))
      stop(sprintf("Covariate '%s' not found in cov.data", col_name))

    col_vals <- cov_data[[col_name]]

    if (is.factor(col_vals)) {
      levs <- levels(col_vals)
      ndum  <- length(levs) - 1L

      if (ndum <= 0) next

      dm <- model.matrix(~ 0 + col_vals)
      colnames(dm) <- levs
      dm <- dm[, levs[-1], drop = FALSE]

      for (j in seq_len(ndum)) {
        coef_name <- paste0(prefix, pidx)
        if (!coef_name %in% names(par))
          stop(sprintf("Parameter '%s' not found for factor '%s' level '%s'",
                       coef_name, col_name, levs[j + 1]))
        result <- result + par[[coef_name]] * dm[, j]
        pidx <- pidx + 1L
      }
    } else {
      coef_name <- paste0(prefix, pidx)
      if (!coef_name %in% names(par))
        stop(sprintf("Parameter '%s' not found for covariate '%s'",
                     coef_name, col_name))
      result <- result + par[[coef_name]] * as.numeric(col_vals)
      pidx <- pidx + 1L
    }
  }

  result
}

pit_linear_fevd <- function(fit, x = NULL, ...) {
  # transform non-linear margins to uniform distribution
  # argument for ...: lower.tail

  if (is.null(x)) x <- fit[["x"]]

  type  <- fit[["type"]]
  par   <- fit[["results"]][["par"]]

  if (type == "GEV") {
    const_loc   <- isTRUE(fit[["const.loc"]])
    const_scale <- isTRUE(fit[["const.scale"]])
    const_shape <- isTRUE(fit[["const.shape"]])

    if (!const_shape) {
      stop("Non-stationary GEV with non-constant shape not yet supported.
  Use pit_gev() with explicit parameter vectors.")
    }
    shape <- par["shape"]

    if (!const_loc) {
      loc <- build_ns_param(
        term_names = fit[["par.models"]][["term.names"]][["location"]],
        cov_data   = fit[["cov.data"]],
        par        = par,
        prefix     = "mu"
      )
    } else {
      loc <- rep(par["location"], length(x))
    }

    if (!const_scale) {
      scale <- build_ns_param(
        term_names = fit[["par.models"]][["term.names"]][["scale"]],
        cov_data   = fit[["cov.data"]],
        par        = par,
        prefix     = "sigma"
      )
    } else {
      scale <- rep(par["scale"], length(x))
    }

    return(pit_gev(x, loc = loc, scale = scale, shape = shape, ...))
  }

  if (type == "GP") {
    idx <- x > fit[["threshold"]]
    x_exc <- x[idx]

    const_scale <- isTRUE(fit[["const.scale"]])
    const_shape <- isTRUE(fit[["const.shape"]])

    if (!const_shape) {
      stop("Non-stationary GP with non-constant shape not yet supported.
  Use pit_gp() with explicit parameter vectors.")
    }
    shape <- par["shape"]

    if (!const_scale) {
      cov_exc <- fit[["cov.data"]][idx, , drop = FALSE]
      scale <- build_ns_param(
        term_names = fit[["par.models"]][["term.names"]][["scale"]],
        cov_data   = cov_exc,
        par        = par,
        prefix     = "sigma"
      )
    } else {
      scale <- rep(par["scale"], length(x_exc))
    }

    return(pit_gp(x_exc, threshold = fit[["threshold"]],
                      scale = scale, shape = shape, ...))
  }

  stop(paste(sprintf("Unknown fevd type: %s", type), "Use 'GEV' or 'GP' only"))
}