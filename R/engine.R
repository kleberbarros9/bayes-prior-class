surv_q_calibration <- function(q) {
  ratio <- -log(.05) / log(2)
  r <- 3.2
  if (q == 1) {
    kappa <- (ratio - r) / (r^2 / (1 + r) - ratio / 2)
    q_one <- 1 + kappa / 2
  } else if (q == 2) {
    kappa <- log(ratio) / log(r)
    q_one <- 1
  } else if (q == 3) {
    kappa <- (ratio - r) / (r^2 - ratio)
    q_one <- 1 + kappa
  } else if (q == 4) {
    kappa <- uniroot(function(k) expm1(k * r) / expm1(k) - ratio,
      c(1e-7, 2), tol = 1e-11)$root
    q_one <- expm1(kappa) / kappa
  } else if (q == 5) {
    kappa <- (ratio - r) / (r * log1p(r) - ratio * log(2))
    q_one <- 1 + kappa * log(2)
  } else {
    return(list(kappa = 1, lambda = 1, baseline_shape = log(r) / log(8)))
  }
  list(kappa = kappa, lambda = log(2) / q_one, baseline_shape = log(r) / log(8))
}


# Select a recognizable calendar unit closest to a maximum of twelve on a ratio scale.
surv_normal_quadrature <- function(points) {
  jacobi <- matrix(0, points, points)
  index <- seq_len(points - 1L)
  jacobi[cbind(index, index + 1L)] <- sqrt(index)
  jacobi[cbind(index + 1L, index)] <- sqrt(index)
  eig <- eigen(jacobi, symmetric = TRUE)
  list(nodes = eig$values, weights = eig$vectors[1, ]^2)
}

surv_stan_data <- function(dataset, model, gh_points) {
  n <- length(dataset$time)
  family <- model$family
  K <- 1L
  M <- I <- exposure <- matrix(0, n, 1)
  reference <- 1
  bin <- rep(1L, n)

  if (family == 5L) {
    cuts <- c(0, .5, 1, 2, 4, Inf)
    K <- length(cuts) - 1L
    exposure <- vapply(seq_len(K), function(k) {
      pmax(0, pmin(dataset$model_time, cuts[k+1]) - cuts[k])
    }, numeric(n))
    reference <- pmax(0, pmin(1, tail(cuts, -1)) - head(cuts, -1))
    bin <- findInterval(dataset$model_time, cuts, left.open = TRUE, all.inside = TRUE)
    M <- I <- matrix(0, n, K)
  }
  if (family == 6L) {
    arguments <- list(knots = c(.5, 1, 2), degree = 3L,
      Boundary.knots = c(0, 20), intercept = TRUE)
    M <- do.call(splines2::mSpline, c(list(x = dataset$model_time), arguments))
    I <- do.call(splines2::iSpline, c(list(x = dataset$model_time), arguments))
    reference <- as.numeric(do.call(splines2::iSpline, c(list(x = 1), arguments)))
    K <- ncol(M)
    exposure <- matrix(0, n, K)
  }
  quadrature <- if (dataset$frailty) surv_normal_quadrature(gh_points) else
    list(nodes = 0, weights = 1)
  calibration <- surv_q_calibration(model$q)

  list(N = n, U = length(dataset$starts), P = ncol(dataset$X),
    J = if (dataset$frailty) length(dataset$starts) else 0L, centered = 0L,
    S = as.integer(family %in% c(2L, 3L, 4L)), K = as.integer(K),
    G = length(quadrature$nodes), family = family, q = model$q,
    time = surv_vector(dataset$time), event = surv_vector(dataset$event, integer = TRUE), X = dataset$X,
    unit_id = surv_vector(dataset$cluster, integer = TRUE),
    unit_start = surv_vector(dataset$starts, integer = TRUE),
    unit_end = surv_vector(dataset$ends, integer = TRUE),
    M = unname(as.matrix(M)), I = unname(as.matrix(I)), exposure = exposure,
    reference = surv_vector(reference), bin = surv_vector(bin, integer = TRUE),
    gh_nodes = surv_vector(quadrature$nodes),
    gh_log_weights = surv_vector(log(quadrature$weights)), kappa = calibration$kappa,
    lambda = calibration$lambda, baseline_shape = calibration$baseline_shape,
    time_ratio = dataset$time_ratio)
}


# Explicit dimensions prevent RStan from treating a one-element vector as a scalar.
surv_vector <- function(x, integer = FALSE) {
  values <- if (integer) as.integer(x) else as.numeric(x)
  array(values, dim = c(length(values)))
}

surv_initial_values <- function(stan_data, chain_id=1L) {
 list(alpha=.02*(chain_id-1), beta=surv_vector(rep(0,stan_data$P)),
 log_shape=surv_vector(rep(0,stan_data$S)),
 increments=surv_vector(rep(0,stan_data$K-1L)),
 smooth_sd=surv_vector(if(stan_data$K>1L) .25 else numeric()))
}
surv_validate_inputs <- function(d) {
 stopifnot(d$N>0, d$P==3, all(is.finite(d$time)), all(d$time>0),
 all(d$event %in% 0:1), identical(dim(d$X),c(as.integer(d$N),as.integer(d$P))),
 nrow(d$M)==d$N, ncol(d$M)==d$K, all(is.finite(d$X)))
 invisible(TRUE)
}

surv_valid_metrics <- function(fit) {
  is.list(fit) && identical(fit$status, "complete") &&
    isTRUE(fit$diagnostics$valid) && is.finite(fit$diagnostics$max_Rhat) &&
    is.finite(fit$diagnostics$min_ESS) &&
    !is.null(fit$pointwise_elpd) && all(is.finite(fit$pointwise_elpd)) &&
    !is.null(fit$loo$estimates) && all(is.finite(fit$loo$estimates)) &&
    !is.null(fit$waic$estimates) && all(is.finite(fit$waic$estimates))
}

surv_parameter_names <- function(stan_data) {
 c("alpha","beta",if(stan_data$S>0) "log_shape",if(stan_data$K>1) c("increments","smooth_sd"))
}
surv_diagnostics <- function(fit, cfg, pars) {
  invalid <- function(reason, sampled = FALSE) list(valid = FALSE, sampled = sampled,
    reason = reason, max_Rhat = NA_real_, min_ESS = NA_real_,
    divergences = NA_integer_, depth_hits = NA_integer_, min_BFMI = NA_real_)
  if (!inherits(fit, "stanfit")) return(invalid("Stan did not return a stanfit object."))
  sampler <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE),
    error = function(error) list())
  if (length(sampler) != cfg$chains || any(!vapply(sampler, function(x) {
    is.matrix(x) && nrow(x) >= 4L && all(c("divergent__", "energy__", "treedepth__") %in% colnames(x))
  }, TRUE))) return(invalid("Stan returned no usable samples for every chain."))
  summary <- tryCatch(rstan::summary(fit, pars = pars)$summary,
    error = function(error) NULL)
  if (is.null(summary) || !nrow(summary)) return(invalid("Empty parameter summary.", TRUE))
  if (any(!is.finite(summary[, c("mean", "sd", "Rhat", "n_eff")]))) {
    return(invalid("Non-finite parameter diagnostics.", TRUE))
  }
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), 0.0))
  depth <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 13), 0.0))
  bfmi <- vapply(sampler, function(x) {
    variance <- var(x[, "energy__"])
    if (!is.finite(variance) || variance <= 0) return(NA_real_)
    mean(diff(x[, "energy__"])^2) / variance
  }, 0.0)
  maximum <- max(summary[, "Rhat"])
  minimum <- min(summary[, "n_eff"])
  valid <- maximum <= 1.01 && minimum >= cfg$min_ess &&
    divergences == 0 && depth == 0 && all(is.finite(bfmi)) && min(bfmi) >= .3
  list(valid = valid, sampled = TRUE, max_Rhat = maximum, min_ESS = minimum,
    divergences = divergences, depth_hits = depth,
    min_BFMI = if (all(is.finite(bfmi))) min(bfmi) else NA_real_,
    parameter_summary = summary)
}

# Valid posterior samples survive a metric failure and can be reused without refitting.
surv_compute_metrics <- function(result, cfg) {
  tryCatch({
    log_lik <- result$log_lik
    dimensions <- dim(log_lik)
    if (length(dimensions) != 3L || any(dimensions < 1L) || any(!is.finite(log_lik))) {
      stop("Pointwise log likelihood must be a finite iteration-by-chain-by-unit array.")
    }
    shifted <- sweep(log_lik, 3L, apply(log_lik, 3L, max), "-")
    relative_eff <- loo::relative_eff(exp(shifted))
    if (length(relative_eff) != dimensions[3] || any(!is.finite(relative_eff)) ||
        any(relative_eff <= 0)) stop("Non-finite relative effective sample sizes.")
    result$loo <- loo::loo(log_lik, r_eff = relative_eff, cores = cfg$cores)
    result$waic <- loo::waic(log_lik)
    result$pareto_k <- loo::pareto_k_values(result$loo)
    result$pointwise_elpd <- result$loo$pointwise[, "elpd_loo"]
    if (any(!is.finite(result$pointwise_elpd)) ||
        any(!is.finite(result$loo$estimates)) || any(!is.finite(result$waic$estimates))) {
      stop("LOO or WAIC returned non-finite values.")
    }
    result$status <- "complete"
    result$error <- NULL
    result
  }, error = function(error) {
    result$status <- "metrics_failed"
    result$error <- conditionMessage(error)
    result
  })
}

surv_fit <- function(compiled, dataset, model, cfg, seed, save_posterior = NULL) {
  stan_data <- surv_stan_data(dataset, model, cfg$gh_points)
  surv_validate_inputs(stan_data)
  history <- list()
  centered <- 0L
  attempts <- cfg$attempts
  for (attempt in seq_len(attempts)) {
    iterations <- cfg$iterations * 2L^(attempt - 1L)
    stan_data$centered <- centered
    cat("  Attempt", attempt, "|", iterations, "iterations per chain",
      "| frailty:", if (stan_data$J == 0) "none" else if (centered) "centered" else "non-centered", "\n")
    flush.console()
    initialize <- function(chain_id = 1L) {
      initial <- surv_initial_values(stan_data, chain_id)
      initial[vapply(initial, length, integer(1)) > 0L]
    }
    fit <- rstan::sampling(compiled, data = stan_data,
      chains = cfg$chains, cores = cfg$cores, iter = iterations,
      warmup = iterations %/% 2L, seed = seed + attempt,
      init = initialize, refresh = max(100L, iterations %/% 4L),
      control = list(adapt_delta = if (attempt == 1L) .99 else .999,
        max_treedepth = 13L))
    diagnostics <- surv_diagnostics(fit, cfg, surv_parameter_names(stan_data))
    diagnostics$parameterization <- if (centered) "centered" else "non-centered"
    history[[attempt]] <- diagnostics
    if (!isTRUE(diagnostics$sampled)) {
      return(list(status = "sampler_failed", fatal = TRUE, diagnostics = diagnostics,
        history = history, error = diagnostics$reason))
    }
    cat(sprintf("  Rhat %.4f | ESS %.0f | divergences %s | depth %s | BFMI %.3f\n",
      diagnostics$max_Rhat, diagnostics$min_ESS, diagnostics$divergences,
      diagnostics$depth_hits, diagnostics$min_BFMI))
    if (isTRUE(diagnostics$valid)) break
    geometry_failed <- !is.finite(diagnostics$min_BFMI) || diagnostics$min_BFMI < .3 ||
      isTRUE(diagnostics$divergences > 0) || isTRUE(diagnostics$depth_hits > 0)
    if (stan_data$J > 0 && geometry_failed) centered <- 1L - centered
    cat("  Retrying with more iterations; parameterization changes only for poor geometry.\n")
  }
  if (!isTRUE(diagnostics$valid)) return(list(status = "diagnostics_failed",
    diagnostics = diagnostics, history = history))

  log_lik <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  dimensions <- dim(log_lik)
  if (length(dimensions) != 3L || dimensions[3] != stan_data$U) {
    stop("Unexpected pointwise log-likelihood dimensions.")
  }
  result <- list(status = "posterior_saved", diagnostics = diagnostics, history = history,
    parameter_summary = rstan::summary(fit, pars = c("alpha", "beta"))$summary,
    log_lik = log_lik, chains = rstan::extract(fit, pars = surv_parameter_names(stan_data), permuted = FALSE), iterations = iterations,
    parameter_draws = as.matrix(fit, pars = surv_parameter_names(stan_data)))
  if (!is.null(save_posterior)) save_posterior(result)
  surv_compute_metrics(result, cfg)
}


