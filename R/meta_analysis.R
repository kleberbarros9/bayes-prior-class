meta_save_rds <- function(object,file,...) {
 tmp<-tempfile("partial_",tmpdir=dirname(file));on.exit(unlink(tmp),add=TRUE)
 base::saveRDS(object,tmp,...)
 if(!file.rename(tmp,file))stop("Checkpoint rename failed: ",file)
}
# Script-relative paths; sourcing loads definitions without starting a long fit.
.univ_root <- local({
  files <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
  if (length(files)) dirname(normalizePath(tail(files, 1L), mustWork = TRUE)) else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE)) else getwd()
  }
})
# Unified Rover analysis: reuse saved metrics and chains; export one interval figure only.
# Existing meta_analysis_results.rds is reused; the active execution call is at the end.

univ_log_sum <- function(x) {
  largest <- max(x)
  if (!is.finite(largest)) return(largest)
  largest + log(sum(exp(x - largest)))
}

univ_progress <- function(label, current, total) {
  cat(sprintf("%s: %d/%d (%.1f%%)\n", label, current, total, 100 * current / total))
  flush.console()
}

# This MAP model follows the binomial-logistic-normal structure in RBesT's AS example.
univ_stan_code <- '
data {
  int<lower=1> J;
  array[J] int<lower=0> y;
  array[J] int<lower=1> n;
  real<lower=0> mu_prior_sd;
  real<lower=0> tau_prior_scale;
  int<lower=0, upper=1> centered;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] z;
}
model {
  mu ~ normal(0, mu_prior_sd);
  tau ~ normal(0, tau_prior_scale);
  if (centered == 1) {
    z ~ normal(mu, tau);
    y ~ binomial_logit(n, z);
  } else {
    z ~ std_normal();
    y ~ binomial_logit(n, mu + tau * z);
  }
}
'

univ_load_data <- function(dataset = "dat.roever2022") {
  original <- read.csv(file.path(Sys.getenv("AJS_PROJECT_ROOT"), "data", "oncology_roever2022.csv"), stringsAsFactors=FALSE)
  original <- original[is.finite(original$dose) & original$dose == 80, , drop = FALSE]
  identity <- paste(original$study, original$year)
  if (anyDuplicated(identity)) stop("Repeated study cohorts at 80 mg/m2 require manual review.")
  data <- data.frame(study = identity, year = original$year,
    y = original$events, n = original$total, dose = original$dose,
    dose_unit = "mg/m2")
  if (nrow(data) != 8L || anyNA(data[c("y", "n")]) ||
      any(data$y < 0 | data$y > data$n | data$n < 1) ||
      any(data$y != floor(data$y) | data$n != floor(data$n))) {
    stop("Unexpected Roever 80 mg/m2 data. Check the installed metadat version.")
  }
  data$y <- as.integer(data$y)
  data$n <- as.integer(data$n)
  rownames(data) <- NULL
  data
}

univ_fit_map <- function(compiled_model, data, cfg, seed) {
  iterations <- cfg$iterations
  warmup <- cfg$warmup
  attempts <- if (cfg$mode == "execution") 3L else 1L
  initial_centered <- as.integer(median(pmin(data$y, data$n - data$y)) >= 20)
  empirical_logit <- qlogis((data$y + .5) / (data$n + 1))
  initial_mu <- median(empirical_logit)
  initial_tau <- max(.15, min(1, sd(empirical_logit)))
  history <- list()
  centered <- initial_centered

  for (attempt in seq_len(attempts)) {
    initialize <- function(chain_id = 1L) {
      mu <- initial_mu + .02 * (chain_id - 1.5)
      tau <- initial_tau * (1 + .05 * chain_id)
      z <- if (centered == 1L) empirical_logit else (empirical_logit - mu) / tau
      list(mu = mu, tau = tau, z = z)
    }
    cat("Historical MAP fit | attempt", attempt, "| iterations", iterations,
      "| parameterization:", if (centered == 1L) "centered" else "non-centered", "\n")
    flush.console()

    fit <- rstan::sampling(
      compiled_model,
      data = list(J = nrow(data), y = data$y, n = data$n,
        mu_prior_sd = cfg$mu_prior_sd, tau_prior_scale = cfg$tau_prior_scale,
        centered = centered),
      chains = cfg$chains,
      iter = iterations,
      warmup = warmup,
      seed = seed + attempt - 1L,
      init = initialize,
      cores = min(cfg$chains, cfg$cores),
      refresh = max(50L, as.integer(iterations / 5L)),
      control = list(adapt_delta = if (attempt == 1L) .99 else .999,
        max_treedepth = 15L)
    )

    summary <- rstan::summary(fit, pars = c("mu", "tau", "z"))$summary
    sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), 0.0))
    depth_hits <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 15), 0.0))
    bfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2) /
      stats::var(x[, "energy__"]), 0.0)
    valid <- all(is.finite(summary[, "Rhat"])) &&
      max(summary[, "Rhat"]) <= 1.01 && min(summary[, "n_eff"]) >= cfg$min_ess &&
      divergences == 0 && depth_hits == 0 && min(bfmi) >= .3

    cat(sprintf(
      "Diagnostics | max Rhat %.5f | min ESS %.1f | divergences %d | depth hits %d | min BFMI %.3f\n",
      max(summary[, "Rhat"]), min(summary[, "n_eff"]), divergences, depth_hits, min(bfmi)
    ))
    flush.console()
    diagnostics <- list(
      parameterization = if (centered == 1L) "centered" else "non-centered",
      valid = valid, parameter_summary = summary, divergences = divergences,
      treedepth_hits = depth_hits, BFMI = bfmi, iterations = iterations,
      warmup = warmup, chains = cfg$chains
    )

    samples <- as.matrix(fit, pars = c("mu", "tau"))
    keep <- unique(round(seq(1, nrow(samples), length.out = min(cfg$map_draws, nrow(samples)))))
    selected <- samples[keep, , drop = FALSE]

    history[[attempt]] <- diagnostics
    if (valid || cfg$mode == "test") {
      diagnostics$attempt_history <- history
      return(list(draws = selected, diagnostics = diagnostics))
    }

    if (attempt < attempts) {
      geometry_failed <- divergences > 0 || depth_hits > 0 ||
        any(!is.finite(bfmi)) || min(bfmi) < .3
      if (geometry_failed) {
        centered <- 1L - centered
        cat("Sampling geometry failed; switching parameterization.\n")
      } else {
        cat("Sampling geometry passed; retaining parameterization and increasing ESS.\n")
      }
      iterations <- 2L * iterations
      warmup <- 2L * warmup
    }
  }

  worst <- rownames(summary)[which.min(summary[, "n_eff"])]
  stop(sprintf(
    paste0("Historical MAP diagnostics failed after %d attempts: max Rhat=%.5f; ",
      "min ESS=%.1f (%s); divergences=%d; depth hits=%d; min BFMI=%.3f. ",
      "Completed blocks remain saved. Send this diagnostic line and warnings()."),
    attempts, max(summary[, "Rhat"]), min(summary[, "n_eff"]), worst,
    divergences, depth_hits, min(bfmi)
  ))
}

# Evaluate the posterior predictive CDF on the log-odds scale in memory-sized blocks.
univ_map_cdf <- function(z, draws) {
  result <- numeric(length(z))
  blocks <- split(seq_along(z), ceiling(seq_along(z) / 256L))
  tau <- pmax(draws[, "tau"], 1e-10)

  for (indices in blocks) {
    standardized <- sweep(outer(z[indices], draws[, "mu"], "-"), 2L, tau, "/")
    result[indices] <- rowMeans(stats::pnorm(standardized))
  }

  result
}

univ_fit_beta <- function(probability_quantiles, probabilities) {
  target <- qlogis(probability_quantiles)
  if (any(!is.finite(target)) || diff(probability_quantiles) <= 0) {
    stop("Invalid target quantiles for Beta calibration.")
  }
  objective <- function(log_parameters) {
    fitted <- qbeta(probabilities, exp(log_parameters[1]), exp(log_parameters[2]))
    if (any(!is.finite(fitted)) || any(fitted <= 0 | fitted >= 1)) return(1e100)
    sum((qlogis(fitted) - target)^2)
  }

  center <- probability_quantiles[1]
  spread <- diff(probability_quantiles)
  approximate_sd <- spread / diff(qnorm(probabilities))
  concentration <- max(.1, center * (1 - center) / approximate_sd^2 - 1)
  calibrated_start <- c(center, 1 - center) * concentration
  starts <- c(lapply(c(.25, 1, 4), function(multiplier) {
    pmax(1e-6, calibrated_start * multiplier)
  }), list(c(1, 3), c(.5, .5), c(2, 5), c(5, 10)))

  fits <- lapply(starts, function(start) optim(
    log(start), objective, method = "L-BFGS-B",
    lower = log(c(1e-8, 1e-8)), upper = log(c(1e14, 1e14)),
    control = list(maxit = 4000, factr = 10, pgtol = 1e-9)
  ))
  best <- fits[[which.min(vapply(fits, function(x) x$value, 0.0))]]
  if (best$value > 1e-7) {
    candidate <- optim(best$par, objective, method = "Nelder-Mead",
      control = list(maxit = 10000, reltol = 1e-12))
    if (candidate$value < best$value) best <- candidate
  }
  parameters <- exp(best$par)
  cdf_error <- max(abs(pbeta(probability_quantiles, parameters[1], parameters[2]) - probabilities))
  if (!is.finite(best$value) || best$value > 1e-7 || cdf_error > 1e-5) {
    stop("Beta two-quantile calibration failed: target = ",
      paste(signif(probability_quantiles, 7), collapse = ", "),
      "; squared logit error = ", best$value, "; CDF error = ", cdf_error)
  }
  parameters
}

# All five maps are increasing from zero to infinity on the nonnegative odds scale.
univ_make_q_map <- function(family, median_log_odds, spread, odds_ratio = 3.2) {
  target_ratio <- -log(.05) / log(2)
  r <- odds_ratio
  baseline_shape <- log(r) / spread

  if (family == "qrog") {
    name <- "QROG"
    parameter <- (target_ratio - r) / (r^2 / (1 + r) - target_ratio / 2)
    log_q <- function(log_u) log_u + log1p(parameter * plogis(log_u))
    formula <- "u + kappa*u^2/(1+u)"
    parameter_name <- "kappa"
  } else if (family == "power") {
    name <- "Power-odds"
    parameter <- log(target_ratio) / log(r)
    log_q <- function(log_u) parameter * log_u
    formula <- "u^kappa"
    parameter_name <- "kappa"
  } else if (family == "quadratic") {
    name <- "Quadratic-odds"
    parameter <- (target_ratio - r) / (r^2 - target_ratio)
    log_q <- function(log_u) log_u + univ_softplus(log(parameter) + log_u)
    formula <- "u + kappa*u^2"
    parameter_name <- "kappa"
  } else if (family == "exponential") {
    name <- "Exponential-odds"

    log_expm1 <- function(x) {
      answer <- numeric(length(x))
      large <- x > 50
      answer[large] <- x[large] + log1p(-exp(-x[large]))
      answer[!large] <- log(expm1(x[!large]))
      answer
    }

    ratio_equation <- function(kappa) {
      log_expm1(kappa * r) - log_expm1(kappa) - log(target_ratio)
    }
    parameter <- uniroot(ratio_equation, c(1e-8, 20), tol = 1e-11)$root
    log_q <- function(log_u) log_expm1(parameter * exp(log_u)) - log(parameter)
    formula <- "expm1(kappa*u)/kappa"
    parameter_name <- "kappa"
  } else if (family == "loglinear") {
    name <- "Loglinear-odds"
    parameter <- (target_ratio - r) /
      (r * log1p(r) - target_ratio * log(2))
    log_q <- function(log_u) log_u + log1p(parameter * univ_softplus(log_u))
    formula <- "u + kappa*u*log1p(u)"
    parameter_name <- "kappa"
  } else {
    stop("Unknown q family: ", family)
  }

  stopifnot(is.finite(parameter), parameter > 0)
  log_q_one <- log_q(0)
  lambda <- exp(log(log(2)) - log_q_one)

  cdf <- function(z) {
    log_u <- baseline_shape * (z - median_log_odds)
    log_hazard <- log(log(2)) + log_q(log_u) - log_q_one
    -expm1(-exp(log_hazard))
  }

  calibration_error <- max(abs(
    cdf(c(median_log_odds, median_log_odds + spread)) - c(.5, .95)
  ))
  if (calibration_error > 1e-7) stop("Two-quantile calibration failed: ", name)

  list(name = name, type = "fixed", cdf = cdf,
    calibration = data.frame(
      model = name, q_formula = formula, shape_name = parameter_name,
      shape_parameter = parameter, lambda = lambda,
      baseline = "Common log-logistic", baseline_shape = baseline_shape,
      baseline_median_odds = exp(median_log_odds), baseline_odds_ratio = r,
      quantile_calibration_error = calibration_error
    ))
}

# Hold G fixed within each fold and calibrate every q map to the same historical quantiles.
univ_models <- function(draws, cfg) {
  map_cdf <- function(z) univ_map_cdf(z, draws)
  quantiles_z <- vapply(c(.5, .95), function(p) {
    uniroot(function(z) map_cdf(z) - p, c(-100, 100), tol = 1e-9)$root
  }, 0.0)
  spread <- diff(quantiles_z)
  stopifnot(spread > 0)
  probabilities <- plogis(quantiles_z)

  proposed <- lapply(cfg$q_families, function(family) {
    univ_make_q_map(family, quantiles_z[1], spread, cfg$baseline_odds_ratio)
  })

  beta_parameters <- univ_fit_beta(probabilities, c(.5, .95))
  beta_cdf <- function(z) pbeta(plogis(z), beta_parameters[1], beta_parameters[2])
  gamma_shape <- exp(uniroot(function(value) {
    shape <- exp(value)
    log(qgamma(.95, shape) / qgamma(.5, shape)) - spread
  }, c(log(.02), log(1e6)), tol = 1e-10)$root)
  gamma_rate <- qgamma(.5, gamma_shape) / exp(quantiles_z[1])
  gamma_cdf <- function(z) pgamma(exp(z), gamma_shape, rate = gamma_rate)
  normal_sigma <- spread / qnorm(.95)
  normal_cdf <- function(z) pnorm(z, quantiles_z[1], normal_sigma)
  robust_cdf <- function(z) cfg$robust_weight * map_cdf(z) +
    (1 - cfg$robust_weight) * plogis(z)

  comparators <- list(
    list(name = "Beta-binomial", type = "fixed", cdf = beta_cdf),
    list(name = "Gamma-odds", type = "fixed", cdf = gamma_cdf),
    list(name = "Lognormal-odds", type = "fixed", cdf = normal_cdf),
    list(name = "MAP", type = "fixed", cdf = map_cdf),
    list(name = "Robust-MAP", type = "fixed", cdf = robust_cdf),
    list(name = "SAM", type = "sam", cdf = robust_cdf)
  )

  models <- c(proposed, comparators)
  names(models) <- vapply(models, function(x) x$name, "")
  calibration <- do.call(rbind, lapply(proposed, function(x) x$calibration))

  list(models = models, calibration = calibration, quantiles_z = quantiles_z,
    quantiles_p = probabilities, beta_parameters = beta_parameters,
    gamma_parameters = c(shape = gamma_shape, rate = gamma_rate),
    normal_parameters = c(mu = quantiles_z[1], sigma = normal_sigma),
    map_cdf = map_cdf)
}

univ_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Exact CDF increments give prior cell masses; mesh refinement controls within-cell error.
univ_prepare_grid <- function(specification, points, cfg) {
  models <- specification$models
  left <- -30
  right <- 30

  for (attempt in seq_len(8L)) {
    omitted <- vapply(models, function(m) m$cdf(left) + 1 - m$cdf(right), 0.0)
    if (max(omitted) < 1e-9) break
    left <- max(-650, left * 1.5)
    right <- min(650, right * 1.5)
  }
  if (max(omitted) >= 1e-9) stop("Prior grid tail bounds were insufficient.")

  edges <- seq(left, right, length.out = points + 1L)
  z <- (head(edges, -1L) + tail(edges, -1L)) / 2
  distributions <- lapply(models, function(m) {
    mass <- pmax(0, diff(m$cdf(edges)))
    mass / sum(mass)
  })
  log_p <- -univ_softplus(-z)
  log_q <- -univ_softplus(z)
  p <- plogis(z)
  center <- sum(distributions$MAP * p)

  list(z = z, p = p, log_p = log_p, log_q = log_q, masses = distributions,
    center = center, omitted = omitted, points = points)
}

univ_sam_weight <- function(y, n, center, delta) {
  alternatives <- c(center - delta, center + delta)
  alternatives <- alternatives[alternatives > 0 & alternatives < 1]
  if (!length(alternatives)) stop("No admissible SAM alternative.")
  log_ratio <- dbinom(y, n, center, log = TRUE) -
    max(dbinom(y, n, alternatives, log = TRUE))
  plogis(log_ratio)
}

univ_prior_mass <- function(name, grid, y, n, cfg) {
  if (name != "SAM") return(grid$masses[[name]])
  weight <- univ_sam_weight(y, n, grid$center, cfg$sam_delta)
  uniform_mass <- diff(plogis(c(
    grid$z[1] - diff(grid$z)[1] / 2,
    grid$z + diff(grid$z)[1] / 2
  )))
  uniform_mass <- uniform_mass / sum(uniform_mass)
  weight * grid$masses$MAP + (1 - weight) * uniform_mass
}

univ_posterior <- function(name, grid, y, n, cfg) {
  prior <- univ_prior_mass(name, grid, y, n, cfg)
  log_likelihood <- y * grid$log_p + (n - y) * grid$log_q
  log_weights <- log(prior) + log_likelihood
  log_normalizer <- univ_log_sum(log_weights)
  weights <- exp(log_weights - log_normalizer)

  list(weights = weights, mean_p = sum(weights * grid$p),
    log_normalizer = log_normalizer, log_likelihood = log_likelihood)
}

univ_count_pmf <- function(weights, grid, n) {
  probability <- vapply(0:n, function(y) {
    sum(weights * dbinom(y, size = n, prob = grid$p))
  }, 0.0)
  probability / sum(probability)
}

univ_predictive_summary <- function(pmf) {
  counts <- seq_along(pmf) - 1L
  cdf <- cumsum(pmf)
  quantile <- function(p) counts[which(cdf >= p)[1]]
  average <- sum(counts * pmf)

  c(mean = average, sd = sqrt(sum((counts - average)^2 * pmf)),
    lower = quantile(.025), median = quantile(.5), upper = quantile(.975))
}

univ_count_summary <- function(weights, grid, n) {
  cdf <- function(k) sum(weights * pbinom(k, n, grid$p))
  quantile <- function(probability) {
    left <- -1
    right <- n
    while (right - left > 1) {
      middle <- floor((left + right) / 2)
      if (cdf(middle) >= probability) right <- middle else left <- middle
    }
    right
  }
  mean_p <- sum(weights * grid$p)
  variance <- sum(weights * n * grid$p * (1 - grid$p)) +
    n^2 * sum(weights * (grid$p - mean_p)^2)
  c(mean = n * mean_p, sd = sqrt(variance), lower = quantile(.025),
    median = quantile(.5), upper = quantile(.975))
}

# Classical criteria use individual Bernoulli outcomes, reconstructed from grouped counts.
univ_criteria <- function(name, grid, y, n, cfg) {
  posterior <- univ_posterior(name, grid, y, n, cfg)
  weights <- posterior$weights
  average <- posterior$mean_p
  variance <- function(x) {
    mean <- sum(weights * x)
    sum(weights * (x - mean)^2)
  }

  lppd <- y * log(average) + (n - y) * log1p(-average)
  p_waic <- y * variance(grid$log_p) + (n - y) * variance(grid$log_q)
  mean_deviance <- -2 * sum(weights * posterior$log_likelihood)
  deviance_mean <- -2 * (y * log(average) + (n - y) * log1p(-average))

  elpd <- 0
  if (y > 0) {
    fit_minus_success <- univ_posterior(name, grid, y - 1L, n - 1L, cfg)
    elpd <- elpd + y * log(fit_minus_success$mean_p)
  }
  if (y < n) {
    fit_minus_failure <- univ_posterior(name, grid, y, n - 1L, cfg)
    elpd <- elpd + (n - y) * log1p(-fit_minus_failure$mean_p)
  }

  fixed_prior <- name != "SAM"
  c(
    ELPD_LOO_patient = elpd,
    LOOIC_patient = -2 * elpd,
    LPML_patient = if (fixed_prior) elpd else NA_real_,
    WAIC_patient = if (fixed_prior) -2 * (lppd - p_waic) else NA_real_,
    p_WAIC = if (fixed_prior) p_waic else NA_real_,
    DIC_patient = if (fixed_prior) 2 * mean_deviance - deviance_mean else NA_real_
  )
}

univ_analyze_fold <- function(specification, target, cfg) {
  previous_signature <- NULL

  for (level in 0:cfg$max_refinements) {
    points <- cfg$grid_points * 2L^level
    cat("  Predictive integration | grid cells:", points, "\n")
    flush.console()
    grid <- univ_prepare_grid(specification, points, cfg)
    summaries <- list()
    pmfs <- list()
    signatures <- list()

    for (index in seq_along(specification$models)) {
      model <- specification$models[[index]]
      name <- model$name
      univ_progress(paste("  Predictive models", name), index, length(specification$models))
      if (is.na(target$y)) {
        if (name == "SAM") {
          pmfs[[name]] <- list(heldout_study = NULL)
          next
        }

        future_pmf <- univ_count_pmf(grid$masses[[name]], grid, target$n)
        pmfs[[name]] <- list(heldout_study = future_pmf)
        signatures[[name]] <- future_pmf
        next
      }

      posterior <- univ_posterior(name, grid, target$y, target$n, cfg)
      predictive <- if (isTRUE(cfg$predictions))
        univ_count_summary(posterior$weights, grid, target$n) else
        setNames(rep(NA_real_, 5L), c("mean", "sd", "lower", "median", "upper"))
      criteria <- univ_criteria(name, grid, target$y, target$n, cfg)

      heldout_score <- if (name != "SAM")
        lchoose(target$n, target$y) +
        univ_log_sum(log(grid$masses[[name]]) + posterior$log_likelihood) else NA_real_
      heldout_interval <- if (isTRUE(cfg$predictions) && name != "SAM") univ_count_summary(grid$masses[[name]], grid, target$n) else rep(NA_real_, 5L)

      summaries[[name]] <- data.frame(
        study = target$study, year = target$year, y = target$y, n = target$n,
        model = name, log_score_heldout_study = heldout_score,
        prior_predictive_mean = heldout_interval[1],
        prior_predictive_lower = heldout_interval[3],
        prior_predictive_upper = heldout_interval[5],
        predictive_replica_mean = predictive["mean"],
        predictive_replica_lower = predictive["lower"],
        predictive_replica_upper = predictive["upper"],
        as.list(criteria), check.names = FALSE
      )
      pmfs[[name]] <- list(heldout_study = NULL, same_study_replica = NULL)
      check_counts <- unique(c(target$y, floor(target$n * c(.025, .1, .25, .5, .75, .9, .975))))
      posterior_cdf <- if (isTRUE(cfg$predictions)) vapply(check_counts, function(k) {
        sum(posterior$weights * pbinom(k, target$n, grid$p))
      }, 0.0) else numeric()
      prior_cdf <- if (isTRUE(cfg$predictions) && name != "SAM") vapply(check_counts, function(k) {
        sum(grid$masses[[name]] * pbinom(k, target$n, grid$p))
      }, 0.0) else numeric()
      signatures[[name]] <- c(
        replica_moments = predictive[c("mean", "sd")][is.finite(predictive[c("mean", "sd")])] / max(1, target$n),
        prior_moments = heldout_interval[1:2][is.finite(heldout_interval[1:2])] / max(1, target$n),
        posterior_cdf = posterior_cdf, prior_cdf = prior_cdf,
        posterior_mean = posterior$mean_p,
        heldout_score = heldout_score[is.finite(heldout_score)],
        criteria = criteria[is.finite(criteria)] / (if (cfg$mode == "test") max(1, target$n) else 1)
      )
    }

    signature <- unlist(signatures, use.names = TRUE)
    changes <- if (is.null(previous_signature)) rep(Inf, length(signature)) else
      abs(signature - previous_signature)
    difference <- max(changes)
    worst_component <- names(signature)[which.max(changes)]
    cat("  Maximum integration change:", format(difference, digits = 5),
      "| component:", worst_component, "\n")
    flush.console()
    if (is.finite(difference) && difference < cfg$grid_tolerance) {
      summary <- if (length(summaries)) do.call(rbind, summaries) else NULL
      if (!is.null(summary)) rownames(summary) <- NULL
      return(list(summary = summary, pmfs = pmfs,
        integration = list(points = points, maximum_difference = difference,
          omitted_prior_mass = grid$omitted),
        historical_quantiles = specification$quantiles_p,
        calibration = specification$calibration,
        conventional_calibration = list(beta = specification$beta_parameters,
          gamma = specification$gamma_parameters,
          lognormal = specification$normal_parameters)))
    }

    previous_signature <- signature
  }

  stop("Predictive grid did not stabilize; maximum difference = ", difference,
    "; component = ", worst_component, "; cells = ", points)
}

# Verify numerical predictive and LOO calculations against an analytical Beta-binomial case.
univ_self_check <- function(cfg) {
  edges <- seq(-30, 30, length.out = 8193L)
  z <- (head(edges, -1L) + tail(edges, -1L)) / 2
  mass <- pmax(0, diff(pbeta(plogis(edges), 2, 5)))
  grid <- list(z = z, p = plogis(z), log_p = -univ_softplus(-z),
    log_q = -univ_softplus(z), masses = list("Beta-binomial" = mass / sum(mass)))
  fitted <- univ_posterior("Beta-binomial", grid, 3L, 12L, cfg)
  computed_pmf <- univ_count_pmf(fitted$weights, grid, 12L)
  k <- 0:12
  exact_pmf <- exp(lchoose(12, k) + lbeta(5 + k, 26 - k) - lbeta(5, 14))
  exact_loo <- 3 * log(4 / 18) + 9 * log(13 / 18)
  computed_loo <- univ_criteria("Beta-binomial", grid, 3L, 12L, cfg)["ELPD_LOO_patient"]
  errors <- c(LOO_error = unname(abs(computed_loo - exact_loo)),
    predictive_PMF_error = max(abs(computed_pmf - exact_pmf)))
  if (max(errors) > 1e-5) stop("Analytical Beta-binomial self-check failed.")
  cat("Analytical Beta-binomial self-check passed.\n")
  errors
}

# Save checkpoints in the same RDS used for final results, keeping three output filenames.
univ_run_dataset <- function(
  mode = c("test", "execution"),
  q_families = c("qrog", "exponential"),
  dataset, output_dir, checkpoint = NULL, historical_rds = NULL, predictions = FALSE
) {
  mode <- match.arg(mode)
  allowed <- c("qrog", "exponential")
  if (!length(q_families) || any(!q_families %in% allowed) || anyDuplicated(q_families)) {
    stop("Choose distinct q families from: ", paste(allowed, collapse = ", "))
  }

  packages <- c("rstan", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing)) {
    stop("Install missing packages first: ", paste(missing, collapse = ", "))
  }

  cfg <- list(
    protocol = "oncology_v1", dataset = dataset, mode = mode, seed = 20261020L,
    q_families = q_families, predictions = predictions, baseline_odds_ratio = 3.2,
    chains = if (mode == "test") 2L else 4L,
    iterations = if (mode == "test") 600L else 4000L,
    warmup = if (mode == "test") 300L else 2000L,
    cores = max(1L, min(4L, as.integer(Sys.getenv("AJS_WORKERS", "4")))),
    map_draws = if (mode == "test") 400L else 8000L,
    grid_points = if (mode == "test") 512L else 2048L,
    max_refinements = 7L, grid_tolerance = if (mode == "test") 5e-5 else 1e-6,
    min_ess = if (mode == "test") 400 else 1000,
    robust_weight = .5, sam_delta = .1, n_future = 30L,
    mu_prior_sd = 10, tau_prior_scale = 1
  )
  numerical_check <- univ_self_check(cfg)
  data <- univ_load_data(dataset)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(output_dir, "results_prior_comparison.rds")
  cat("Mode:", mode, "\nOutput directory:", normalizePath(output_dir), "\n")
  print(data)
  flush.console()

  notes <- c(
    "Roever 2022: eight study cohorts at the fixed dose of 80 mg/m2; retrospective leave-one-study-out validation.",
    "This dataset, dose and q choices follow exploratory comparisons; this is higher-precision reanalysis, not independent confirmation.",
    "Dose schedules and clinical populations may differ; a common dose does not remove all clinical heterogeneity.",
    "QROG and Exponential-odds share a fixed log-logistic baseline within each fold and historical MAP Q50/Q95.",
    "Beta, Gamma-odds and Lognormal-odds match those historical quantiles; MAP mixtures are not quantile matched.",
    "Hyperpriors: mu ~ Normal(0,10^2), tau ~ HalfNormal(0,1). Full mode uses four chains and up to 8000 retained MAP draws.",
    "Study scores evaluate a new study; patient criteria evaluate a new patient within a study and condition on historical information.",
    "LPML equals exact patient LOO ELPD for fixed priors; LOOIC=-2*ELPD. These are not independent measures of evidence.",
    "SAM is data-adaptive and refitted in patient LOO; study scores, WAIC, DIC and LPML are omitted for SAM.",
    "Paired patient ELPD uncertainty is clustered by study; eight overlapping validation folds limit the SE interpretation.",
    paste("Optional predictions:", if (predictions) "enabled" else "disabled; fit criteria retained"),
    "Test mode is only a workflow check. Full-mode numerical convergence uses unscaled fit criteria."
  )
  references <- c(
    "Neuenschwander et al. (2010), Clinical Trials 7:5-18. DOI:10.1177/1740774509356002.",
    "Ursino et al. (2021): DOI:10.1214/20-AOAS1390.",
    "Roever et al. (2022): DOI:10.1002/sim.9484.",
    "Schmidli et al. (2014), Biometrics 70:1023-1032. DOI:10.1111/biom.12242.",
    "Weber et al. (2021), Journal of Statistical Software 100(19). DOI:10.18637/jss.v100.i19.",
    "Yang et al. (2023), Biometrics 79:2857-2868. DOI:10.1111/biom.13927.",
    "Vehtari, Gelman and Gabry (2017), Statistics and Computing 27:1413-1432. DOI:10.1007/s11222-016-9696-4."
  )
  results <- list(status = "partial", configuration = cfg, data = data, folds = list(),
    numerical_check = numerical_check,
    notes = notes, references = references, sessionInfo = sessionInfo())

  if (file.exists(rds_path)) {
    previous <- readRDS(rds_path)
    if (is.null(previous$configuration$predictions)) previous$configuration$predictions <- TRUE
    compatible_fields <- setdiff(names(cfg), c("max_refinements", "predictions"))
    if (identical(previous$configuration[compatible_fields], cfg[compatible_fields]) &&
        identical(previous$data, data)) {
      results <- previous
      results$configuration <- cfg
      cat("Compatible checkpoint found; completed fits will be reused.\n")
    }
  }

  started <- proc.time()[3]
  fold_indices <- if (mode == "test") which.min(data$n) else seq_len(nrow(data))
  required <- c(as.character(fold_indices), if (predictions) "future")
  source_folds <- list()
  if (!is.null(historical_rds) && file.exists(historical_rds)) {
    source <- readRDS(historical_rds)
    compatible <- identical(source$data, data) &&
      identical(source$configuration$mu_prior_sd, cfg$mu_prior_sd) &&
      identical(source$configuration$tau_prior_scale, cfg$tau_prior_scale)

    if (compatible) {
      for (key in intersect(required, names(source$folds))) {
        saved <- source$folds[[key]]
        expected_studies <- if (key == "future") data$study else data$study[-as.integer(key)]
        usable <- identical(saved$historical_studies, expected_studies) &&
          !is.null(saved$MAP_draws) && isTRUE(saved$diagnostics$valid)
        if (usable) source_folds[[key]] <- saved
      }
      cat("Reusable historical MAP blocks:", length(source_folds), "\n")
    }
  }

  missing_blocks <- setdiff(required, names(results$folds))
  needs_fit <- any(!missing_blocks %in% names(source_folds))
  compiled <- NULL

  if (needs_fit) {
    cat("Compiling the historical MAP model; the first compilation can take a few minutes.\n")
    flush.console()
    compiled <- rstan::stan_model(model_code = univ_stan_code, model_name = "historical_MAP_binomial")
  }

  for (position in seq_along(required)) {
    key <- required[position]
    univ_progress("Analysis blocks", position, length(required))

    if (key %in% names(results$folds)) {
      cat("  Reusing saved block:", key, "\n")
      next
    }

    if (key == "future") {
      historical <- data
      target <- data.frame(study = "Future study", year = NA_real_, y = NA_integer_, n = cfg$n_future)
    } else {
      index <- as.integer(key)
      historical <- data[-index, , drop = FALSE]
      target <- data[index, , drop = FALSE]
      cat("  Held-out study:", target$study, "| y =", target$y, "| n =", target$n, "\n")
    }

    if (key %in% names(source_folds)) {
      saved <- source_folds[[key]]
      fit <- list(draws = saved$MAP_draws, diagnostics = saved$diagnostics)
      cat("  Reusing historical MAP draws; Stan is not repeated.\n")
    } else {
      fit <- univ_fit_map(compiled, historical, cfg, cfg$seed + position)
    }

    if (nrow(fit$draws) > cfg$map_draws) {
      selected <- unique(round(seq(1, nrow(fit$draws), length.out = cfg$map_draws)))
      fit$draws <- fit$draws[selected, , drop = FALSE]
    }
    specification <- univ_models(fit$draws, cfg)
    analysis <- univ_analyze_fold(specification, target, cfg)

    if (key == "future") {
      analysis$summary <- NULL
      analysis$pmfs <- lapply(analysis$pmfs, function(x) list(heldout_study = x$heldout_study))
    }

    results$folds[[key]] <- c(analysis, list(
      MAP_draws = fit$draws, diagnostics = fit$diagnostics,
      historical_studies = historical$study, target = target
    ))
    results$status <- "partial"
    meta_save_rds(results, rds_path, compress = FALSE)
    if (!is.null(checkpoint)) checkpoint(results)
    cat("  Checkpoint saved.\n")
    flush.console()
  }

  results$status <- if (mode == "test") "test_complete" else "complete"
  results$last_run_seconds <- unname(proc.time()[3] - started)
  results <- univ_finalize(results)
  if (!is.null(checkpoint)) checkpoint(results)
  invisible(results)
}

univ_finalize <- function(results) {
  folds <- results$folds[setdiff(names(results$folds), "future")]
  summaries <- do.call(rbind, lapply(folds, function(x) x$summary))
  rownames(summaries) <- NULL
  results$study_results <- summaries
  models <- unique(summaries$model)
  baseline <- summaries[summaries$model == "Beta-binomial", ]

  comparison <- lapply(models, function(name) {
    model <- summaries[summaries$model == name, ]
    model <- model[match(baseline$study, model$study), ]
    differences <- model$log_score_heldout_study - baseline$log_score_heldout_study
    valid <- is.finite(differences)
    total <- if (all(valid)) sum(model$log_score_heldout_study) else NA_real_
    paired_se <- if (sum(valid) > 1) sqrt(sum(valid) * var(differences[valid])) else NA_real_
    data.frame(
      model = name, study_folds = nrow(model), ELPD_heldout_study = total,
      delta_ELPD_vs_Beta = if (all(valid)) sum(differences) else NA_real_,
      SE_delta_ELPD = paired_se,
      ELPD_LOO_patient = sum(model$ELPD_LOO_patient),
      LOOIC_patient = sum(model$LOOIC_patient),
      LPML_patient = if (all(is.finite(model$LPML_patient))) sum(model$LPML_patient) else NA_real_,
      WAIC_patient = if (all(is.finite(model$WAIC_patient))) sum(model$WAIC_patient) else NA_real_,
      DIC_patient = if (all(is.finite(model$DIC_patient))) sum(model$DIC_patient) else NA_real_,
      predictive_count_coverage = if (all(valid)) mean(
        model$y >= model$prior_predictive_lower & model$y <= model$prior_predictive_upper
      ) else NA_real_
    )
  })
  results$comparison <- do.call(rbind, comparison)
  candidate_names <- results$folds[[names(folds)[1]]]$calibration$model
  q_ranking <- results$comparison[results$comparison$model %in% candidate_names, ]
  q_ranking <- q_ranking[order(-q_ranking$ELPD_heldout_study), ]
  reference <- summaries[summaries$model == "QROG", ]

  q_ranking$delta_ELPD_vs_QROG <- NA_real_
  q_ranking$SE_delta_vs_QROG <- NA_real_
  if (nrow(reference)) {
    for (index in seq_len(nrow(q_ranking))) {
      candidate <- summaries[summaries$model == q_ranking$model[index], ]
      candidate <- candidate[match(reference$study, candidate$study), ]
      difference <- candidate$log_score_heldout_study - reference$log_score_heldout_study
      q_ranking$delta_ELPD_vs_QROG[index] <- sum(difference)
      q_ranking$SE_delta_vs_QROG[index] <- if (length(difference) > 1L)
        sqrt(length(difference) * var(difference)) else NA_real_
    }
  }
  results$q_ranking <- q_ranking
  references <- intersect(c("QROG", "Exponential-odds"), models)
  competing <- setdiff(models, references)
  rows <- list()
  for (proposed in references) {
    left <- summaries[summaries$model == proposed, ]
    for (competitor in competing) {
      right <- summaries[summaries$model == competitor, ]
      right <- right[match(left$study, right$study), ]
      study_difference <- left$log_score_heldout_study - right$log_score_heldout_study
      patient_difference <- left$ELPD_LOO_patient - right$ELPD_LOO_patient
      paired_se <- function(x) {
        if (any(!is.finite(x)) || length(x) < 2L) return(NA_real_)
        sqrt(length(x) * var(x))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        proposed = proposed, comparator = competitor,
        delta_study_ELPD = if (all(is.finite(study_difference))) sum(study_difference) else NA_real_,
        SE_delta_study = paired_se(study_difference),
        delta_patient_ELPD = sum(patient_difference),
        SE_delta_patient_by_study = paired_se(patient_difference),
        delta_WAIC = sum(left$WAIC_patient) - sum(right$WAIC_patient),
        delta_DIC = sum(left$DIC_patient) - sum(right$DIC_patient)
      )
    }
  }
  results$proposed_comparisons <- do.call(rbind, rows)
  results
}

# Re-export saved results without refitting historical models or repeating validation.
univ_export_dataset <- function(results, output_dir = NULL) {
  if (is.character(results)) {
    path <- results
    results <- readRDS(path)
    if (is.null(output_dir)) output_dir <- dirname(path)
  }
  if (!results$status %in% c("complete", "test_complete")) {
    stop("Checkpoint is incomplete; call run_univariate_rover() to resume.")
  }
  if (is.null(output_dir)) {
    output_dir <- file.path(file.path(.univ_root, "results"), "meta_analysis", "oncology")
    if (results$configuration$mode == "test") output_dir <- file.path(output_dir, "test")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  meta_save_rds(results, file.path(output_dir, "results_prior_comparison.rds"), compress = FALSE)

  report <- capture.output({
    cat("META-ANALYSIS: PRIOR COMPARISON\nStatus:", results$status, "\n")
    print(results$configuration)
    print(results$numerical_check)
    cat("\nDATA\n")
    print(results$data, row.names = FALSE)
    cat("\nNOTES\n")
    cat(results$notes, sep = "\n")
    cat("\n\nCROSS-VALIDATED COMPARISON\n")
    print(results$comparison, row.names = FALSE, digits = 7)
    cat("\nQ-MAP RANKING (EXPLORATORY)\n")
    print(results$q_ranking, row.names = FALSE, digits = 7)
    cat("\nPAIRED COMPARISONS: POSITIVE ELPD AND NEGATIVE WAIC/DIC FAVOR PROPOSED\n")
    print(results$proposed_comparisons, row.names = FALSE, digits = 7)
    cat("\nSTUDY-SPECIFIC RESULTS\n")
    print(results$study_results, row.names = FALSE, digits = 7)

    for (key in names(results$folds)) {
      fold <- results$folds[[key]]
      cat("\nBLOCK:", key, "\n")
      print(fold$target)
      cat("Historical predictive Q50/Q95 on probability scale:\n")
      print(fold$historical_quantiles)
      print(fold$calibration, row.names = FALSE)
      print(fold$conventional_calibration)
      print(fold$integration)
      print(fold$diagnostics)
    }

    cat("\nREFERENCES\n")
    cat(results$references, sep = "\n")
    cat("\n\nLast run seconds:", results$last_run_seconds, "\n")
    print(results$sessionInfo)
  })
  writeLines(report, file.path(output_dir, "report_prior_comparison.txt"), useBytes = TRUE)
  cat("\nSaved numerical results and report in:", normalizePath(output_dir), "\n")
  invisible(results)
}

# The active call at the end runs the full metrics analysis.
# run_univariate_rover("test") or run_univariate_rover("execution")


run_rover_analysis <- function(
  mode = c("test", "exec"),
  prediction = NULL,
  datasets = "dat.roever2022",
  q_families = c("qrog", "exponential"),
  historical_rds = NULL
) {
  mode <- match.arg(mode, c("test", "exec", "execution"))
  if (mode == "exec") mode <- "execution"
  if (!is.null(prediction) && !identical(prediction, "pred")) {
    stop("The second argument must be omitted or equal to 'pred'.")
  }
  predictions <- identical(prediction, "pred")
  allowed <- "dat.roever2022"
  if (!length(datasets) || any(!datasets %in% allowed) || anyDuplicated(datasets)) {
    stop("Select distinct datasets from: ", paste(allowed, collapse = ", "))
  }
  output_dir <- file.path(file.path(.univ_root, "results"), "meta_analysis", "oncology")
  if (mode == "test") output_dir <- file.path(output_dir, "test")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(output_dir, "meta_analysis_results.rds")
  bundle <- list(protocol = "oncology_v1", mode = mode,
    q_families = q_families, predictions = predictions, datasets = list(), status = "partial")

  if (file.exists(rds_path)) {
    previous <- readRDS(rds_path)
    if (identical(previous$protocol, bundle$protocol) &&
        identical(previous$mode, mode) && identical(previous$q_families, q_families)) {
      bundle <- previous
    }
  }
  bundle$predictions <- predictions
  cat("Optional predictions:", if (predictions) "enabled" else "disabled (fast mode)", "\n")
  historical_bundle <- NULL
  if (!is.null(historical_rds) && file.exists(historical_rds)) {
    historical_bundle <- readRDS(historical_rds)
    cat("Historical bundle available; valid MAP fits will be reused.\n")
  }
  work_dir <- tempfile("meta_analysis_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  for (index in seq_along(datasets)) {
    dataset <- datasets[index]
    cat("\nDATASET", index, "OF", length(datasets), ":", dataset, "\n")
    current <- univ_load_data(dataset)
    saved <- bundle$datasets[[dataset]]
    if (!is.null(saved) && !identical(saved$data, current)) saved <- NULL
    has_predictions <- !is.null(saved$folds$future) &&
      !identical(saved$configuration$predictions, FALSE)
    if (!is.null(saved) && saved$status %in% c("complete", "test_complete") &&
        (!predictions || has_predictions)) {
      cat("Reusing completed dataset.\n")
      next
    }
    internal_path <- file.path(work_dir, "results_prior_comparison.rds")
    if (file.exists(internal_path)) unlink(internal_path)
    if (!is.null(saved)) meta_save_rds(saved, internal_path)

    checkpoint <- function(result) {
      bundle$datasets[[dataset]] <<- result
      bundle$status <<- "partial"
      meta_save_rds(bundle, rds_path, compress = FALSE)
    }
    historical_path <- NULL
    old <- if (!is.null(saved)) saved else historical_bundle$datasets[[dataset]]
    if (!is.null(old)) {
      historical_path <- file.path(work_dir, "historical_source.rds")
      meta_save_rds(old, historical_path)
    }
    result <- univ_run_dataset(mode, q_families, dataset, work_dir,
      checkpoint, historical_rds = historical_path, predictions = predictions)
    bundle$datasets[[dataset]] <- result
    meta_save_rds(bundle, rds_path, compress = FALSE)
  }
  bundle$requested_datasets <- datasets
  bundle$status <- if (mode == "test") "test_complete" else "complete"
  meta_save_rds(bundle, rds_path, compress = FALSE)
  export_univariate_rover(bundle, output_dir)
  invisible(bundle)
}

export_univariate_rover <- function(results, output_dir = NULL) {
  if (is.character(results)) {
    if (is.null(output_dir)) output_dir <- dirname(results)
    results <- readRDS(results)
  }
  if (is.null(output_dir)) {
    output_dir <- file.path(file.path(.univ_root, "results"), "meta_analysis", "oncology")
    if (results$mode == "test") output_dir <- file.path(output_dir, "test")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  meta_save_rds(results, file.path(output_dir, "meta_analysis_results.rds"), compress = FALSE)
  keys <- results$requested_datasets
  if (is.null(keys)) keys <- names(results$datasets)
  completed <- keys[vapply(results$datasets[keys], function(x) {
    !is.null(x) && x$status %in% c("complete", "test_complete")
  }, TRUE)]
  if (!length(completed)) stop("No completed dataset available to export.")
  report <- c("ROEVER 2022: 80 mg/m2 - TWO GENERATED PRIORS AND SIX COMPARATORS",
    paste("Status:", results$status),
    "All selected datasets must be reported, including unfavorable results.")
  work_dir <- tempfile("meta_report_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  for (key in completed) {
    item <- results$datasets[[key]]
    if (identical(results$predictions, FALSE)) item$configuration$predictions <- FALSE
    univ_export_dataset(item, work_dir)
    report <- c(report, "", paste("DATASET:", key),
      readLines(file.path(work_dir, "report_prior_comparison.txt"), warn = FALSE))
  }
  writeLines(report, file.path(output_dir, "meta_analysis_report.txt"), useBytes = TRUE)
  cat("\nCombined outputs:", normalizePath(output_dir), "\n")
  invisible(results)
}




# A deterministic normal expectation preserves the original SAM weight without a dense grid.
univ_map_center <- function(draws) {
  previous <- NULL
  for (points in c(32L, 64L, 128L, 256L, 512L)) {
    jacobi <- matrix(0, points, points)
    indices <- seq_len(points - 1L)
    jacobi[cbind(indices, indices + 1L)] <- sqrt(indices)
    jacobi[cbind(indices + 1L, indices)] <- sqrt(indices)
    normal <- eigen(jacobi, symmetric = TRUE)
    values <- outer(normal$values, pmax(draws[, "tau"], 1e-10))
    values <- sweep(values, 2L, draws[, "mu"], "+")
    center <- sum(normal$vectors[1, ]^2 * rowMeans(plogis(values)))
    if (!is.null(previous) && abs(center - previous) < 1e-8) return(center)
    previous <- center
  }
  stop("SAM historical mean did not stabilize; saved results remain unchanged.")
}

univ_prior_random <- function(size, name, fold, cfg, center) {
  conventional <- fold$conventional_calibration
  map_random <- function(size) {
    selected <- sample.int(nrow(fold$MAP_draws), size, replace = TRUE)
    draws <- fold$MAP_draws[selected, , drop = FALSE]
    plogis(rnorm(size, draws[, "mu"], pmax(draws[, "tau"], 1e-10)))
  }

  if (name %in% c("QROG", "Exponential-odds")) {
    calibration <- fold$calibration[fold$calibration$model == name, ]
    if (nrow(calibration) != 1L) stop("Missing calibration: ", name)
    kappa <- calibration$shape_parameter
    exponential <- rexp(size, rate = calibration$lambda)
    if (name == "QROG") {
      root <- sqrt((1 - exponential)^2 + 4 * (1 + kappa) * exponential)
      u <- 2 * exponential / (root + 1 - exponential)
    } else {
      u <- log1p(kappa * exponential) / kappa
    }
    return(plogis(log(calibration$baseline_median_odds) +
      log(u) / calibration$baseline_shape))
  }

  if (name == "Beta-binomial") {
    return(rbeta(size, conventional$beta[1], conventional$beta[2]))
  }
  if (name == "Gamma-odds") {
    odds <- rgamma(size, shape = conventional$gamma["shape"],
      rate = conventional$gamma["rate"])
    return(plogis(log(odds)))
  }
  if (name == "Lognormal-odds") {
    return(plogis(rnorm(size, conventional$lognormal["mu"],
      conventional$lognormal["sigma"])))
  }
  if (name == "MAP") return(map_random(size))
  if (name %in% c("Robust-MAP", "SAM")) {
    weight <- if (name == "SAM") univ_sam_weight(
      fold$target$y, fold$target$n, center, cfg$sam_delta
    ) else cfg$robust_weight
    from_map <- runif(size) < weight
    values <- runif(size)
    values[from_map] <- map_random(sum(from_map))
    return(values)
  }
  stop("Unsupported model: ", name)
}

# Rejection sampling uses the exact binomial likelihood envelope at y/n.
univ_posterior_random <- function(name, fold, cfg, center, size, seed) {
  set.seed(seed)
  y <- fold$target$y
  n <- fold$target$n
  log_bound <- dbinom(y, n, y / n, log = TRUE)
  accepted <- numeric()
  proposals <- 0L
  total_accepted <- 0L

  while (length(accepted) < size) {
    batch_size <- max(4096L, 2L * (size - length(accepted)))
    candidates <- univ_prior_random(batch_size, name, fold, cfg, center)
    log_acceptance <- dbinom(y, n, candidates, log = TRUE) - log_bound
    keep <- log(runif(batch_size)) <= pmin(0, log_acceptance)
    proposals <- proposals + batch_size
    total_accepted <- total_accepted + sum(keep)
    accepted <- c(accepted, candidates[keep])
    cat(sprintf("    %-18s %d/%d draws | proposals %d | acceptance %.1f%%\n",
      name, min(size, length(accepted)), size, proposals, 100 * total_accepted / proposals))
    flush.console()
    if (proposals > max(2000000, 200 * size) && length(accepted) < size) {
      stop("Low rejection acceptance for ", name, "; completed blocks remain saved.")
    }
  }

  samples <- accepted[seq_len(size)]
  quantiles <- quantile(samples, c(.025, .5, .975), names = FALSE)
  blocks <- split(samples, rep(seq_len(10L), length.out = size))
  block_quantiles <- vapply(blocks, function(x) {
    quantile(x, c(.025, .5, .975), names = FALSE)
  }, numeric(3))
  mcse <- apply(block_quantiles, 1L, sd) / sqrt(length(blocks))
  list(draws = samples, lower = quantiles[1], estimate = quantiles[2],
    upper = quantiles[3], mcse_lower = mcse[1], mcse_median = mcse[2],
    mcse_upper = mcse[3], proposals = proposals,
    acceptance = total_accepted / proposals)
}


# Shared hyperparameters induce posterior dependence between study-specific probabilities.
univ_joint_exponential <- function(data, cfg) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan for the pooled joint fit.")
  calibration <- univ_make_q_map("exponential", 0, qnorm(.95),
    cfg$baseline_odds_ratio)$calibration
  kappa <- calibration$shape_parameter
  lambda <- calibration$lambda
  code <- '
  data {
    int<lower=1> J;
    array[J] int<lower=0> y;
    array[J] int<lower=1> n;
    real<lower=0> kappa;
    real<lower=0> lambda;
    real<lower=0> factor_scale;
    real<lower=0> mu_sd;
    real<lower=0> tau_scale;
  }
  parameters {
    real mu;
    real<lower=0> tau;
    vector[J] log_e;
  }
  transformed parameters {
    vector[J] eta;
    for (j in 1:J) {
      real log_u = log(log1p_exp(log(kappa) + log_e[j])) - log(kappa);
      eta[j] = mu + tau * factor_scale * log_u;
    }
  }
  model {
    mu ~ normal(0, mu_sd);
    tau ~ normal(0, tau_scale);
    target += sum(log_e) - lambda * sum(exp(log_e));
    y ~ binomial_logit(n, eta);
  }
  generated quantities {
    vector[J] p = inv_logit(eta);
    real pooled = 0;
    real total_n = 0;
    for (j in 1:J) {
      pooled += n[j] * p[j];
      total_n += n[j];
    }
    pooled /= total_n;
  }
  '
  cat("\nFitting joint exponential-q hierarchy for the pooled probability.\n")
  compiled <- rstan::stan_model(model_code = code, model_name = "joint_exponential_q_rover")
  iterations <- if (cfg$mode == "test") 1000L else 4000L
  chains <- if (cfg$mode == "test") 2L else 4L
  history <- list()

  for (attempt in seq_len(3L)) {
    cat("Joint pooled fit | attempt", attempt, "| iterations", iterations, "\n")
    fit <- rstan::sampling(compiled,
      data = list(J = nrow(data), y = data$y, n = data$n,
        kappa = kappa, lambda = lambda,
        factor_scale = qnorm(.95) / log(cfg$baseline_odds_ratio),
        mu_sd = cfg$mu_prior_sd, tau_scale = cfg$tau_prior_scale),
      chains = chains, iter = iterations, warmup = iterations %/% 2L,
      cores = min(chains, cfg$cores), seed = cfg$seed + 90000L + attempt,
      init = function() list(mu = qlogis((sum(data$y) + .5) / (sum(data$n) + 1)),
        tau = .5, log_e = rep(log(log(2) / lambda), nrow(data))),
      refresh = max(100L, iterations %/% 5L),
      control = list(adapt_delta = .999, max_treedepth = 15L))
    summary <- rstan::summary(fit, pars = c("mu", "tau", "log_e", "p", "pooled"))$summary
    sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), 0.0))
    depth <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 15), 0.0))
    bfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2) / var(x[, "energy__"]), 0.0)
    min_ess <- if (cfg$mode == "test") 100 else 1000
    valid <- all(is.finite(summary[, "Rhat"])) && max(summary[, "Rhat"]) <= 1.01 &&
      min(summary[, "n_eff"]) >= min_ess && divergences == 0 && depth == 0 &&
      all(is.finite(bfmi)) && min(bfmi) >= .3
    diagnostics <- list(valid = valid, summary = summary, divergences = divergences,
      depth_hits = depth, BFMI = bfmi, iterations = iterations, chains = chains)
    history[[attempt]] <- diagnostics
    cat(sprintf("Joint diagnostics | Rhat %.5f | ESS %.1f | divergences %d | depth %d | BFMI %.3f\n",
      max(summary[, "Rhat"]), min(summary[, "n_eff"]), divergences, depth, min(bfmi)))
    if (valid) {
      samples <- as.matrix(fit, pars = c("mu", "tau", "p", "pooled"))
      interval <- quantile(samples[, "pooled"], c(.025, .5, .975), names = FALSE)
      return(list(draws = samples,
        chains = as.array(fit, pars = c("mu", "tau", "p", "pooled")),
        diagnostics = diagnostics, history = history,
        lower = interval[1], estimate = interval[2], upper = interval[3],
        calibration = calibration, weights = data$n / sum(data$n),
        definition = "sum(n_j * p_j) / sum(n_j), from a joint exponential-q hierarchy",
        model_note = "Additional joint inferential model; not the leave-one-study-out model used for ranking."))
    }
    iterations <- iterations * 2L
  }
  stop("Joint pooled diagnostics failed. Previous study fits and metrics remain saved.")
}

univ_export_mc_intervals <- function(results, output_dir) {
  data <- results$data
  total_y <- sum(data$y)
  total_n <- sum(data$n)
  exact <- binom.test(total_y, total_n)
  empirical <- data.frame(
    study = "All studies pooled", model = "Observed pooled proportion",
    y = total_y, n = total_n, lower = unname(exact$conf.int[1]),
    estimate = total_y / total_n, upper = unname(exact$conf.int[2]),
    interval = "Exact binomial 95% CI", grid_points = NA_integer_, maximum_change = NA_real_
  )
  keys <- setdiff(names(results$folds), "future")
  intervals <- lapply(keys, function(key) {
    fold <- results$folds[[key]]
    do.call(rbind, lapply(names(fold$posterior_mc$samples), function(name) {
      value <- fold$posterior_mc$samples[[name]]
      data.frame(study = fold$target$study, model = name, y = fold$target$y,
        n = fold$target$n, lower = value$lower, estimate = value$estimate,
        upper = value$upper, interval = "Monte Carlo 95% CrI",
        grid_points = NA_integer_, maximum_change = NA_real_)
    }))
  })
  bayesian <- do.call(rbind, intervals)
  order_models <- unique(bayesian$model)
  bayesian$model <- factor(bayesian$model, levels = rev(order_models))
  model_colors <- c(
    "QROG" = "#0072B2", "Exponential-odds" = "#D55E00",
    "Beta-binomial" = "#009E73", "Gamma-odds" = "#CC79A7",
    "Lognormal-odds" = "#8C6D00", "MAP" = "#000000",
    "Robust-MAP" = "#56B4E9", "SAM" = "#7B3294"
  )
  model_lines <- c(
    "QROG" = "solid", "Exponential-odds" = "longdash",
    "Beta-binomial" = "dashed", "Gamma-odds" = "dotdash",
    "Lognormal-odds" = "twodash", "MAP" = "dotted",
    "Robust-MAP" = "4313", "SAM" = "F282"
  )
  study_order <- unique(bayesian$study)
  titles <- vapply(study_order, function(study) {
    row <- bayesian[bayesian$study == study, ][1, ]
    paste0(study, " (", row$y, "/", row$n, ")")
  }, "")
  bayesian$panel <- factor(bayesian$study, levels = study_order, labels = titles)

  joint <- results$joint_exponential
  joint_row <- empirical
  joint_row$model <- "Exponential q: joint pooled posterior"
  joint_row$lower <- joint$lower
  joint_row$estimate <- joint$estimate
  joint_row$upper <- joint$upper
  joint_row$interval <- "Joint posterior equal-tailed 95% CrI"
  pooled_rows <- rbind(empirical, joint_row)
  pooled_rows$position <- c(2, 1)
  pooled_rows$group <- c("Binomial", "Exponential q")
  pooled_rows$label <- sprintf("%.3f [%.3f, %.3f]",
    pooled_rows$estimate, pooled_rows$lower, pooled_rows$upper)

  top <- ggplot2::ggplot(pooled_rows,
    ggplot2::aes(x = estimate, y = position, colour = group, linetype = group)) +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = position), linewidth = 1.3) +
    ggplot2::geom_point(ggplot2::aes(shape = group), size = 3) +
    ggplot2::geom_text(ggplot2::aes(x = .98, label = label),
      hjust = 1, colour = "black", size = 3.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("Binomial" = "black", "Exponential q" = "#D55E00")) +
    ggplot2::scale_linetype_manual(values = c("Binomial" = "solid", "Exponential q" = "longdash")) +
    ggplot2::scale_shape_manual(values = c("Binomial" = 15, "Exponential q" = 16)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::scale_y_continuous(breaks = c(1, 2), limits = c(.4, 2.6),
      labels = c("Joint exponential q", "Observed binomial")) +
    ggplot2::labs(title = "Pooled DLT probability: binomial reference and joint Bayesian estimate",
      subtitle = "Exact binomial 95% CI versus sample-size-weighted posterior 95% CrI",
      x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11))

  bottom <- ggplot2::ggplot(bayesian,
    ggplot2::aes(x = estimate, y = model, colour = model, linetype = model)) +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = model),
      linewidth = 1) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::facet_wrap(~panel, ncol = 2) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::scale_colour_manual(values = model_colors,
      breaks = order_models, name = NULL) +
    ggplot2::scale_linetype_manual(values = model_lines,
      breaks = order_models, name = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "Study-specific posterior medians and 95% credible intervals",
      x = "Probability of dose-limiting toxicity", y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(1.5, "cm"),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )

  pdf_path <- file.path(output_dir, "meta_analysis_intervals.pdf")
  csv_path <- file.path(output_dir, "posterior_intervals.csv")
  table <- rbind(pooled_rows[names(empirical)],
    transform(bayesian[names(empirical)], model = as.character(model)))
  utils::write.csv(table, csv_path, row.names = FALSE)

  grDevices::pdf(pdf_path, width = 10, height = 13, useDingbats = FALSE)
  tryCatch({
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = 2, ncol = 1, heights = grid::unit(c(.16, .84), "npc")
    )))
    print(top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(bottom, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::popViewport()
  }, finally = grDevices::dev.off())

  cat("\nFigure:", pdf_path, "\nIntervals:", csv_path, "\n")
  invisible(list(empirical = empirical, bayesian = bayesian))
}


# Completed metrics and MAP chains are reused; only absent Monte Carlo blocks are generated.
run_univariate_rover <- function(
  mode = "exec", prediction = NULL,
  datasets = "dat.roever2022", q_families = c("qrog", "exponential"),
  historical_rds = NULL, mc_draws = NULL, refresh_mc = FALSE
) {
  mode <- match.arg(mode, c("exec", "test", "execution"))
  if (mode == "execution") mode <- "exec"
  if (!identical(datasets, "dat.roever2022") ||
      !identical(q_families, c("qrog", "exponential"))) {
    stop("This final script uses Rover and the two prespecified q functions.")
  }
  if (!is.null(prediction) && !identical(prediction, "pred")) {
    stop("The second argument must be omitted or 'pred'.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  if (is.null(mc_draws)) mc_draws <- if (mode == "test") 2000L else 40000L
  if (length(mc_draws) != 1L || !is.finite(mc_draws) || mc_draws < 1000 ||
      mc_draws != floor(mc_draws)) stop("mc_draws must be an integer >= 1000.")

  output_dir <- file.path(
    .univ_root,
    "results/meta_analysis/oncology"
  )
  if (mode == "test") output_dir <- file.path(output_dir, "test")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Keep the existing checkpoint name to reuse previous completed analyses.
  path <- file.path(output_dir, "meta_analysis_results.rds")
  bundle <- if (file.exists(path)) readRDS(path) else NULL
  item <- bundle$datasets[[datasets]]
  completed <- !is.null(item) && item$status %in% c("complete", "test_complete")
  expected_mode <- if (mode == "exec") "execution" else "test"
  if (completed && (!identical(item$configuration$mode, expected_mode) ||
      !identical(item$configuration$q_families, q_families))) {
    stop("Saved configuration differs. Keep test and execution results in separate folders.")
  }
  if (completed && any(!vapply(item$folds, function(fold) {
    !is.null(fold$MAP_draws) && isTRUE(fold$diagnostics$valid)
  }, TRUE))) stop("Saved MAP diagnostics are not valid; inspect the checkpoint.")

  need_prediction <- identical(prediction, "pred") && is.null(item$folds$future)
  if (!completed || need_prediction) {
    cat("Resuming missing analysis blocks; saved historical chains will be reused.\n")
    bundle <- run_rover_analysis(mode, prediction, datasets, q_families, historical_rds)
    item <- bundle$datasets[[datasets]]
  } else {
    cat("SKIP: completed dataset, historical chains and evaluation metrics.\n")
    export_univariate_rover(bundle, output_dir)
  }

  keys <- setdiff(names(item$folds), "future")
  models <- c("QROG", "Exponential-odds", "Beta-binomial", "Gamma-odds",
    "Lognormal-odds", "MAP", "Robust-MAP", "SAM")
  for (position in seq_along(keys)) {
    key <- keys[position]
    fold <- item$folds[[key]]
    signature <- list(version = "posterior_rejection_v1", draws = mc_draws,
      MAP_draws = fold$MAP_draws, target = fold$target,
      calibration = fold$calibration, conventional = fold$conventional_calibration,
      robust_weight = item$configuration$robust_weight, sam_delta = item$configuration$sam_delta)
    cached <- fold$posterior_mc
    if (!refresh_mc && identical(cached$signature, signature) &&
        all(models %in% names(cached$samples))) {
      cat("SKIP Monte Carlo:", fold$target$study, "\n")
      next
    }
    if (is.null(cached) || !identical(cached$signature, signature) || refresh_mc) {
      cached <- list(signature = signature, samples = list())
    }
    center <- univ_map_center(fold$MAP_draws)
    cat("\nStudy", position, "of", length(keys), ":", fold$target$study, "\n")
    for (index in seq_along(models)) {
      name <- models[index]
      if (!is.null(cached$samples[[name]])) {
        cat("  SKIP saved samples:", name, "\n")
        next
      }
      cached$samples[[name]] <- univ_posterior_random(name, fold,
        item$configuration, center, mc_draws,
        item$configuration$seed + 1000L * as.integer(key) + index)
      item$folds[[key]]$posterior_mc <- cached
      bundle$datasets[[datasets]] <- item
      temporary <- tempfile("checkpoint_", tmpdir = output_dir)
      meta_save_rds(bundle, temporary, compress = FALSE)
      if (!file.rename(temporary, path)) stop("Could not save Monte Carlo checkpoint.")
    }
  }
  joint_signature <- list(version = "joint_exponential_q_v1", data = item$data,
    mode = item$configuration$mode, mu_sd = item$configuration$mu_prior_sd,
    tau_scale = item$configuration$tau_prior_scale,
    baseline_odds_ratio = item$configuration$baseline_odds_ratio)
  if (!identical(item$joint_exponential$signature, joint_signature) ||
      !isTRUE(item$joint_exponential$diagnostics$valid)) {
    item$joint_exponential <- univ_joint_exponential(item$data, item$configuration)
    item$joint_exponential$signature <- joint_signature
    bundle$datasets[[datasets]] <- item
    temporary <- tempfile("checkpoint_", tmpdir = output_dir)
    meta_save_rds(bundle, temporary, compress = FALSE)
    if (!file.rename(temporary, path)) stop("Could not save joint pooled fit.")
  } else {
    cat("SKIP: valid cached joint exponential-q pooled fit.\n")
  }
  univ_export_mc_intervals(item, output_dir)
  diagnostics <- do.call(rbind, lapply(keys, function(key) {
    fold <- item$folds[[key]]
    do.call(rbind, lapply(models, function(name) {
      value <- fold$posterior_mc$samples[[name]]
      data.frame(study = fold$target$study, model = name,
        draws = length(value$draws), acceptance = value$acceptance,
        mcse_lower = value$mcse_lower, mcse_median = value$mcse_median,
        mcse_upper = value$mcse_upper)
    }))
  }))
  report_path <- file.path(output_dir, "meta_analysis_report.txt")
  old_report <- readLines(report_path, warn = FALSE)
  marker <- which(old_report == "MONTE CARLO POSTERIOR INTERVALS")
  if (length(marker)) old_report <- head(old_report, marker[1] - 1L)
  writeLines(c(old_report, "MONTE CARLO POSTERIOR INTERVALS",
    "Posterior rejection samples conditional on the saved historical MAP draws.",
    "MCSE estimates quantify simulation error conditional on the saved MAP distribution.",
    capture.output(print(diagnostics, row.names = FALSE, digits = 5)),
    "JOINT EXPONENTIAL-Q POOLED ESTIMATE",
    item$joint_exponential$definition, item$joint_exponential$model_note,
    capture.output(print(item$joint_exponential[c("lower", "estimate", "upper", "diagnostics")]))), report_path)
  cat("\nDone. Saved metric values are unchanged.\n")
  invisible(bundle)
}

