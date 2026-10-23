survival_log_sum <- function(values) {
  largest <- max(values)
  if (!is.finite(largest)) return(largest)
  largest + log(sum(exp(values - largest)))
}

survival_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

survival_splines <- function(time) {
  arguments <- list(knots = c(.5, 1, 2), degree = 3L,
    Boundary.knots = c(0, 20), intercept = TRUE)
  clamped <- pmax(0, pmin(20, time))
  M <- do.call(splines2::mSpline, c(list(x = clamped), arguments))
  I <- do.call(splines2::iSpline, c(list(x = clamped), arguments))
  M_end <- as.numeric(do.call(splines2::mSpline, c(list(x = 20), arguments)))
  I <- I + outer(pmax(time - 20, 0), M_end)
  list(M = as.matrix(M), I = as.matrix(I),
    reference = as.numeric(do.call(splines2::iSpline, c(list(x = 1), arguments))))
}

survival_profile <- function(parameters, family) {
  shape <- if ("log_shape[1]" %in% names(parameters)) exp(parameters["log_shape[1]"]) else 1
  if (family %in% c(5L, 6L)) {
    increments <- parameters[grep("^increments\\[", names(parameters))]
    profile <- c(0, cumsum(parameters["smooth_sd[1]"] * increments))
    weights <- if (family == 6L) exp(profile - survival_log_sum(profile)) else exp(profile)
    reference <- if (family == 6L) survival_splines(1)$reference else c(.5, .5, 0, 0, 0)
    normalization <- sum(weights * reference)
  } else {
    weights <- 1
    normalization <- 1
  }
  list(shape = as.numeric(shape), weights = as.numeric(weights), normalization = as.numeric(normalization))
}

survival_baseline <- function(time, profile, family, basis = NULL) {
  positive <- time > 0
  H <- rep(0, length(time))
  log_h <- rep(NA_real_, length(time))
  t <- time[positive]
  if (!length(t)) return(list(H = H, log_h = log_h))
  shape <- profile$shape
  if (family %in% c(1L, 2L)) {
    H[positive] <- t^shape
    log_h[positive] <- log(shape) + (shape - 1) * log(t)
  } else if (family == 5L) {
    cuts <- c(0, .5, 1, 2, 4, Inf)
    exposure <- vapply(seq_along(profile$weights), function(k) {
      pmax(0, pmin(t, cuts[k + 1L]) - cuts[k])
    }, numeric(length(t)))
    exposure <- matrix(exposure, nrow = length(t))
    bin <- findInterval(t, cuts, left.open = TRUE, all.inside = TRUE)
    H[positive] <- as.numeric(exposure %*% profile$weights) / profile$normalization
    log_h[positive] <- log(profile$weights[bin]) - log(profile$normalization)
  } else if (family == 6L) {
    if (is.null(basis)) basis <- survival_splines(time)
    H <- as.numeric(basis$I %*% profile$weights) / profile$normalization
    log_h[positive] <- log(as.numeric(basis$M[positive, , drop = FALSE] %*% profile$weights)) -
      log(profile$normalization)
  }

  list(H = H, log_h = log_h)
}

survival_survival <- function(parameters, model, dataset, time, basis = NULL) {
  profile <- survival_profile(parameters, model$family)
  eta <- as.numeric(dataset$X %*% parameters[paste0("beta[", seq_len(ncol(dataset$X)), "]")])
  reference_time <- time / dataset$reference_years
  if (model$family == 3L) {
    z <- outer(reference_time, -parameters["alpha"] - eta,
      function(t, mu) (log(t) - mu) / profile$shape)
    return(matrix(pnorm(z, lower.tail = FALSE), nrow = length(time)))
  }
  if (model$family == 4L) {
    z <- profile$shape * outer(log(reference_time), parameters["alpha"] + eta, "+")
    return(exp(-survival_softplus(z)))
  }
  baseline <- survival_baseline(reference_time, profile, model$family, basis)
  exp(-outer(baseline$H, exp(parameters["alpha"] + eta)))
}

survival_loglik <- function(parameters, model, dataset) {
  profile <- survival_profile(parameters, model$family)
  eta <- as.numeric(dataset$X %*% parameters[paste0("beta[", seq_len(ncol(dataset$X)), "]")])
  t <- dataset$model_time
  event <- dataset$event
  if (model$family == 3L) {
    mu <- -parameters["alpha"] - eta
    log_f <- dnorm(log(t), mu, profile$shape, log = TRUE) - log(t)
    log_s <- pnorm(log(t), mu, profile$shape, lower.tail = FALSE, log.p = TRUE)
    return(ifelse(event == 1, log_f, log_s) + event * log(dataset$time_ratio))
  }
  if (model$family == 4L) {
    z <- profile$shape * (log(t) + parameters["alpha"] + eta)
    log_f <- log(profile$shape) - log(t) + z - 2 * survival_softplus(z)
    return(ifelse(event == 1, log_f, -survival_softplus(z)) + event * log(dataset$time_ratio))
  }
  baseline <- survival_baseline(t, profile, model$family)
  log_h <- parameters["alpha"] + eta + baseline$log_h
  H <- exp(parameters["alpha"] + eta) * baseline$H
  ifelse(event == 1, log_h, 0) - H + event * log(dataset$time_ratio)
}


survival_inverse <- function(target, profile, family, spline_grid = NULL) {
  if (family %in% c(1L, 2L)) return(target^(1 / profile$shape))
  if (family == 5L) {
    cuts <- c(0, .5, 1, 2, 4)
    lengths <- diff(c(cuts, Inf))
    cumulative <- c(0, cumsum(profile$weights[1:4] * lengths[1:4]))
    raw_target <- target * profile$normalization
    bin <- findInterval(raw_target, cumulative, all.inside = FALSE)
    bin <- pmax(1L, pmin(length(cuts), bin))
    return(cuts[bin] + (raw_target - cumulative[bin]) / profile$weights[bin])
  }
  if (family == 6L) {
    H <- as.numeric(spline_grid$basis$I %*% profile$weights) / profile$normalization
    end_hazard <- sum(spline_grid$basis$M[nrow(spline_grid$basis$M), ] * profile$weights) /
      profile$normalization
    inverse <- approx(H, spline_grid$time, xout = pmin(target, tail(H, 1)), ties = "ordered")$y
    tail <- target > tail(H, 1)
    inverse[tail] <- 20 + (target[tail] - tail(H, 1)) / end_hazard
    check <- survival_baseline(inverse, profile, family)$H
    error <- max(abs(check - target) / pmax(1, target))
    if (!is.finite(error)) stop("Non-finite spline inversion check.")
    if (error > 1e-4) {
      inverse <- vapply(target, function(value) {
        if (value > tail(H, 1)) return(20 + (value - tail(H, 1)) / end_hazard)
        uniroot(function(t) survival_baseline(t, profile, family)$H - value,
          c(0, 20), tol = 1e-9)$root
      }, 0.0)
    }
    return(inverse)
  }
  stop("Unsupported inverse baseline.")
}
