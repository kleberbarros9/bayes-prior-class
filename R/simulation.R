# General QROG validation. Call run_qrog("test") or run_qrog("execution").
# Outputs: results/simulation/ with section-specific English filenames.

run_qrog <- function(mode = c("test", "execution")) {
  mode <- match.arg(mode)
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop('Install first: install.packages("ggplot2")')
  # Test mode checks the workflow; execution mode runs the full experiment.
  cfg <- list(
    M = 5000L, M_censoring = 1000L, R_sbc = 5000L, L_sbc = 100L,
    grid_points = 2049L, grid_tolerance = 2e-4,
    exposures = c(5, 20, 100, 500), true_rates = c(.05, .2, .4, .8),
    median = .2, upper = .8, probabilities = c(.5, .95),
    threshold = .4, horizon = 5, future_exposure = 10,
    mixture_weight = .5, delta_sam = .1,
    baseline_shapes = c(1, .85), sdlog_baseline = 2,
    seed = 20261019L, N_rejection = 5000L,
    mode = mode,
    censoring_sizes = c(20, 100, 500),
    censoring_fractions = c(.2, .8),
    censoring_rates = c(.2, .4),
    benchmark_rates = c(.05, .2, .8)
  )
  if (mode == "test") {
    cfg$M <- 5L
    cfg$M_censoring <- 3L
    cfg$R_sbc <- 3L
    cfg$L_sbc <- 10L
    cfg$N_rejection <- 20L
    cfg$exposures <- 5
    cfg$true_rates <- .2
    cfg$censoring_sizes <- 10
    cfg$censoring_fractions <- .8
    cfg$censoring_rates <- .2
    cfg$benchmark_rates <- .2
  }
  
  output_dir <- file.path(getwd(), "results", "simulation")
  if (mode == "test") {
    output_dir <- file.path(output_dir, "test")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Mode:", mode, "\n")
  cat("Output directory:", normalizePath(output_dir), "\n")
  flush.console()
  start_time <- proc.time()[3]
  RNGkind("L'Ecuyer-CMRG")
  set.seed(cfg$seed)
  warnings_log <- character()
  old_options <- options(warn = 1)
  on.exit(options(old_options), add = TRUE)
  
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  dir.create(checkpoint_dir, showWarnings = FALSE)
  checkpoint_key <- digest::digest(list(cfg, readLines(file.path(Sys.getenv("AJS_PROJECT_ROOT"), "R", "simulation.R"))), algo="sha256")
  stage_read <- function(name) {
    path <- file.path(checkpoint_dir, paste0(name, ".rds"))
    if (!file.exists(path)) return(NULL)
    value <- readRDS(path)
    if (!identical(value$key, checkpoint_key)) stop("Simulation checkpoint configuration changed: ", path)
    message("Reusing simulation stage: ", name)
    value$values
  }
  stage_save <- function(name, variables, frame) {
    variables <- unique(c(variables, "n_refinements", "n_adaptive_fallbacks", "max_difference", "warnings_log"))
    values <- mget(variables, envir=frame)
    temporary <- tempfile(tmpdir=checkpoint_dir)
    saveRDS(list(key=checkpoint_key, values=values), temporary, compress=FALSE)
    if (!file.rename(temporary, file.path(checkpoint_dir, paste0(name, ".rds")))) stop("Checkpoint rename failed")
  }
  log_warning <- function(x) {
    warnings_log <<- c(warnings_log, x)
    message(x)
  }
  
  bind_rows <- function(x) {
    if (length(x)) do.call(rbind, x) else data.frame()
  }
  
  mcse <- function(x) {
    if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  }
  
  lse <- function(x) {
    a <- max(x)
    if (!is.finite(a)) return(a)
    a+log(sum(exp(x-a)))
  }
  
  log_add <- function(a, b) {
    m <- pmax(a, b)
    m+log(exp(a-m)+exp(b-m))
  }
  
  trap <- function(x, y) sum(diff(x)*(head(y, -1)+tail(y, -1))/2)
  
  make_progress <- function(label, total, every = 100L) {
    started <- proc.time()[3]
    
    function(current) {
      if (current %in% c(0L, 1L, total) || current %% every == 0L) {
        elapsed <- unname(proc.time()[3] - started)
        remaining <- if (current > 0) elapsed * (total - current) / current else NA_real_
        
        cat(sprintf(
          "  %s: %d/%d (%.1f%%) | elapsed %.1fs | ETA %s\n",
          label, current, total, 100 * current / max(total, 1L), elapsed,
          if (is.finite(remaining)) sprintf("%.1fs", remaining) else "--"
        ))
        
        flush.console()
      }
    }
  }
  
  qinv <- function(u, k) {
    v <- numeric(length(u))
    small <- u < 1
    v[small] <- 2*u[small]/(1-u[small]+
                              sqrt((1-u[small])^2+4*(1+k)*u[small]))
    v[!small] <- (u[!small]-1+
                    sqrt((1-u[!small])^2+4*(1+k)*u[!small]))/(2*(1+k))
    v
  }
  
  qfun <- function(u, k) u+k*u*(u/(1+u))
  
  qp <- function(u, k) 1+k-k/(1+u)^2
  # Compute baseline odds and their Jacobian without subtracting probabilities.
  
  base_ll <- function(beta, scale) list(
    name = paste0("LL", beta), scale = scale,
    logodds = function(t) beta*(log(t)-log(scale)),
    logjac = function(t) log(beta)-log(scale)+(beta-1)*(log(t)-log(scale)),
    inverse = function(u) scale*exp(log(u)/beta))
  
  base_ln <- function(s, scale) list(
    name = paste0("LN", s), scale = scale,
    logodds = function(t) {
      z <- (log(t)-log(scale))/s
      pnorm(z, log.p=TRUE)-pnorm(z, lower.tail=FALSE, log.p=TRUE)
    },
    logjac = function(t) {
      z <- (log(t)-log(scale))/s
      dnorm(z, log=TRUE)-log(t)-log(s)-
        2*pnorm(z, lower.tail=FALSE, log.p=TRUE)
    },
    inverse = function(u) {
      z <- log(u)-log1p(u)
      scale*exp(s*qnorm(z, log.p=TRUE))
    })
  
  calibrate_qrog <- function(base) {
    r <- exp(base$logodds(c(cfg$median, cfg$upper)))
    A <- r*r/(1+r)
    constants <- -log1p(-cfg$probabilities)
    R <- constants[2]/constants[1]
    if (!(r[2]/r[1] <= R && R < A[2]/A[1]))
      stop("Inadmissible quantiles for baseline ", base$name)
    k <- (R*r[1]-r[2])/(A[2]-R*A[1])
    lam <- constants[1]/qfun(r[1], k)
    list(name=paste0("QROG-", base$name), type="numerical", base=base,
         lambda=lam, kappa=k,
         logp=function(t) {
           u <- exp(base$logodds(t))
           ans <- log(lam)+log(qp(u, k))+base$logjac(t)-lam*qfun(u, k)
           ans[!is.finite(u)] <- -Inf
           ans
         },
         quantile=function(p) base$inverse(qinv(-log1p(-p)/lam, k)),
         sample=function(n) base$inverse(qinv(rexp(n, lam), k)))
  }
  bases <- c(lapply(cfg$baseline_shapes, base_ll, scale=cfg$median),
             list(base_ln(cfg$sdlog_baseline, cfg$median)))
  proposed_models <- lapply(bases, calibrate_qrog)
  
  fit_gamma <- function(x, p) {
    target_ratio <- x[2]/x[1]
    
    f <- function(z) {
      a <- exp(z)
      log(qgamma(p[2], a))-log(qgamma(p[1], a))-log(target_ratio)
    }
    a <- exp(uniroot(f, c(log(.02), log(10000)), tol=1e-11)$root)
    c(shape=a, rate=qgamma(p[1], a)/x[1])
  }
  # Match Gamma, Lognormal and Weibull to the same two elicited quantiles.
  gg <- unname(fit_gamma(c(cfg$median, cfg$upper), cfg$probabilities))
  mu <- log(cfg$median)
  sig <- (log(cfg$upper)-mu)/qnorm(cfg$probabilities[2])
  w_shape <- log(log1p(-cfg$probabilities[2])/log1p(-cfg$probabilities[1]))/
    log(cfg$upper/cfg$median)
  w_scale <- cfg$median/(-log1p(-cfg$probabilities[1]))^(1/w_shape)
  
  gamma_model <- function(name, a, b) list(name=name, type="gamma", a=a, b=b,
                                           logp=function(t) dgamma(t, a, rate=b, log=TRUE),
                                           quantile=function(p) qgamma(p, a, rate=b), sample=function(n) rgamma(n, a, rate=b))
  gam <- gamma_model("Gamma-Poisson", gg[1], gg[2])
  ln <- list(name="Lognormal-Poisson", type="numerical",
             logp=function(t) dlnorm(t, mu, sig, log=TRUE),
             quantile=function(p) qlnorm(p, mu, sig), sample=function(n) rlnorm(n, mu, sig))
  wei <- list(name="Weibull-Poisson", type="numerical",
              logp=function(t) dweibull(t, w_shape, w_scale, log=TRUE),
              quantile=function(p) qweibull(p, w_shape, w_scale),
              sample=function(n) rweibull(n, w_shape, w_scale))
  # Fixed robustification settings; no historical meta-analysis is fitted.
  a0 <- .25
  b0 <- a0/(gg[1]/gg[2])
  mixture <- list(name="Robust-mixture", type="mixture", w=cfg$mixture_weight,
                  a=c(gg[1], a0), b=c(gg[2], b0))
  # SAM uses the published likelihood-ratio rule, followed by mixture updating.
  sam <- mixture
  sam$name <- "SAM-Poisson"
  sam$type <- "sam"
  sam$theta_h <- gg[1]/gg[2]
  sam$delta <- cfg$delta_sam
  stopifnot(sam$theta_h > sam$delta)
  
  mixture$logp <- function(t) log_add(log(cfg$mixture_weight)+gam$logp(t),
                                      log1p(-cfg$mixture_weight)+dgamma(t, a0, rate=b0, log=TRUE))
  models <- c(proposed_models, list(gam, ln, wei, mixture, sam))
  names(models) <- vapply(models, `[[`, "", "name")
  fixed_models <- models[vapply(models, function(m)m$type!="sam", TRUE)]
  # Integrate on log(theta), including the Jacobian and adaptive tail bounds.
  
  grid <- function(m, a, b, N=cfg$grid_points) {
    
    lf <- function(z) {
      t <- exp(z)
      val <- a*z-b*t+m$logp(t)+z
      val[!is.finite(val)] <- -Inf
      val
    }
    guess <- log(max((a+.5)/(b+1), cfg$median/100))
    center <- log(cfg$median)
    lo <- min(center-35, guess-20)
    hi <- max(center+12, guess+12)
    opt <- optimize(function(z)-lf(z), c(lo, hi))
    mode <- opt$minimum
    peak <- lf(mode)
    if (!is.finite(peak)) stop("Non-finite mode: ", m$name)
    left <- right <- mode
    # Expand geometrically to cover slowly decaying tails on the log scale.
    step <- .25
    
    for (j in seq_len(100L)) {
      left <- max(-700, left-step)
      if (lf(left)<peak-42) break
      if (left<=-700) stop("Lower tail exceeds the numerical range: ", m$name)
      step <- step*1.4
    }
    if (lf(left)>=peak-42) stop("Lower bound not found: ", m$name)
    step <- .25
    
    for (j in seq_len(100L)) {
      right <- min(700, right+step)
      if (lf(right)<peak-42) break
      if (right>=700) stop("Upper tail exceeds the numerical range: ", m$name)
      step <- step*1.4
    }
    if (lf(right)>=peak-42) stop("Upper bound not found: ", m$name)
    z <- seq(left, right, length.out=N)
    z <- sort(unique(c(z, log(cfg$threshold))))
    z <- z[z>=left & z<=right]
    t <- exp(z)
    v <- exp(lf(z)-peak)
    area <- diff(z)*(head(v, -1)+tail(v, -1))/2
    Z <- sum(area)
    cdf <- c(0, cumsum(area))/Z
    v <- v/Z
    
    q <- function(p) exp(approx(cdf, z, xout=p, ties="ordered", rule=2)$y)
    
    f <- function(x) approx(z, cdf, xout=log(x), rule=2)$y
    mean <- trap(z, t*v)
    m2 <- trap(z, t*t*v)
    pred <- trap(z, exp(-cfg$horizon*t)*v)
    pred0 <- trap(z, exp(-cfg$future_exposure*t)*v)
    c(mean=mean, sd=sqrt(max(0, m2-mean^2)), lower=q(.025),
      median=q(.5), upper=q(.975), tail=1-f(cfg$threshold),
      survival=pred, prob_zero=pred0, logZ=log(Z)+peak) -> summary
    list(summary=summary, z=z, density_z=v, cdf=cdf, q=q, F=f, lf=lf, peak=peak)
  }
  n_refinements <- 0L
  n_adaptive_fallbacks <- 0L
  max_difference <- 0
  
  adaptive_posterior <- function(ob) {
    lower_bound <- min(ob$z)
    upper_bound <- max(ob$z)
    
    kernel <- function(z) exp(ob$lf(z) - ob$peak)
    
    integrate_kernel <- function(fun, lower = lower_bound, upper = upper_bound) {
      if (upper <= lower) return(0)
      
      integrate(
        fun, lower, upper,
        rel.tol = 1e-9, abs.tol = 1e-12, subdivisions = 2000L
      )$value
    }
    
    normalizer <- integrate_kernel(kernel)
    
    expectation <- function(fun) {
      integrate_kernel(function(z) fun(exp(z)) * kernel(z)) / normalizer
    }
    
    posterior_cdf <- function(x) {
      vapply(x, function(value) {
        if (value <= exp(lower_bound)) return(0)
        if (value >= exp(upper_bound)) return(1)
        
        integrate_kernel(kernel, upper = log(value)) / normalizer
      }, numeric(1))
    }
    
    posterior_quantile <- function(p) {
      vapply(p, function(probability) {
        exp(uniroot(
          function(z) posterior_cdf(exp(z)) - probability,
          c(lower_bound, upper_bound), tol = 1e-10
        )$root)
      }, numeric(1))
    }
    
    posterior_mean <- expectation(identity)
    posterior_second_moment <- expectation(function(x) x^2)
    
    ob$summary <- c(
      mean = posterior_mean,
      sd = sqrt(max(0, posterior_second_moment - posterior_mean^2)),
      lower = posterior_quantile(.025),
      median = posterior_quantile(.5),
      upper = posterior_quantile(.975),
      tail = 1 - posterior_cdf(cfg$threshold),
      survival = expectation(function(x) exp(-cfg$horizon * x)),
      prob_zero = expectation(function(x) exp(-cfg$future_exposure * x)),
      logZ = log(normalizer) + ob$peak
    )
    
    ob$q <- posterior_quantile
    ob$F <- posterior_cdf
    ob$expectation <- expectation
    ob$method <- "adaptive_quadrature"
    ob
  }
  
  numerical <- function(m, a, b, return_grid = FALSE) {
    previous <- grid(m, a, b)
    
    # Refine without relaxing tolerance; use adaptive quadrature if needed.
    for (level in seq_len(5L)) {
      points <- 2L^level * (cfg$grid_points - 1L) + 1L
      current <- grid(m, a, b, N = points)
      scale <- pmax(abs(current$summary), .01)
      difference <- max(abs(previous$summary - current$summary) / scale)
      max_difference <<- max(max_difference, difference)
      
      if (is.finite(difference) && difference <= cfg$grid_tolerance) {
        current$method <- "refined_grid"
        if (return_grid) return(current)
        return(current$summary)
      }
      
      n_refinements <<- n_refinements + 1L
      previous <- current
    }
    
    n_adaptive_fallbacks <<- n_adaptive_fallbacks + 1L
    cat("  Adaptive quadrature:", m$name, "| a =", a, "| b =", b, "\n")
    flush.console()
    current <- adaptive_posterior(previous)
    
    if (return_grid) current else current$summary
  }
  
  posterior_gamma <- function(a, b) {
    c(mean=a/b, sd=sqrt(a)/b, lower=qgamma(.025, a, rate=b),
      median=qgamma(.5, a, rate=b), upper=qgamma(.975, a, rate=b),
      tail=pgamma(cfg$threshold, a, rate=b, lower.tail=FALSE),
      survival=exp(a*(log(b)-log(b+cfg$horizon))),
      prob_zero=exp(a*(log(b)-log(b+cfg$future_exposure))), logZ=NA_real_)
  }
  
  posterior_mix <- function(m, a, b) {
    w <- m$w
    if(m$type=="sam") {
      
      ll <- function(t) a*log(t)-b*t
      logR <- ll(m$theta_h)-max(ll(m$theta_h-m$delta), ll(m$theta_h+m$delta))
      w <- plogis(logR)
    }
    lm <- lgamma(m$a+a)-lgamma(m$a)+m$a*log(m$b)-(m$a+a)*log(m$b+b)
    lw <- c(log(w), log1p(-w))+lm
    pp <- exp(lw-lse(lw))
    ap <- m$a+a
    bp <- m$b+b
    
    F <- function(t) sum(pp*pgamma(t, ap, rate=bp))
    
    Q <- function(p) {
      hi <- max(qgamma(p, ap, rate=bp))
      uniroot(function(t) F(t)-p, c(0, hi*1.001), tol=1e-10)$root
    }
    mean <- sum(pp*ap/bp)
    c(mean=mean, sd=sqrt(max(0, sum(pp*ap*(ap+1)/bp^2)-mean^2)),
      lower=Q(.025), median=Q(.5), upper=Q(.975),
      tail=sum(pp*pgamma(cfg$threshold, ap, rate=bp, lower.tail=FALSE)),
      survival=sum(pp*exp(ap*(log(bp)-log(bp+cfg$horizon)))),
      prob_zero=sum(pp*exp(ap*(log(bp)-log(bp+cfg$future_exposure)))),
      logZ=lse(lw), prior_weight=w, posterior_weight=pp[1])
  }
  
  posterior <- function(m, a, b) {
    if(m$type=="gamma") posterior_gamma(m$a+a, m$b+b)
    else if(m$type %in% c("mixture", "sam")) posterior_mix(m, a, b)
    else numerical(m, a, b)
  }
  
  references <- c(
    "Yang P, Zhao Y, Nie L, Vallejo J, Yuan Y (2023). SAM: Self-adapting mixture prior to dynamically borrow information from historical data in clinical trials. Biometrics 79:2857-2868. DOI:10.1111/biom.13927. https://pmc.ncbi.nlm.nih.gov/articles/PMC10842647/",
    "Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust meta-analytic-predictive priors in clinical trials with historical control information. Biometrics 70:1023-1032. DOI:10.1111/biom.12242.",
    "Weber S, Li Y, Seaman JW III, Kakizume T, Schmidli H (2021). Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools. Journal of Statistical Software 100(19):1-32. DOI:10.18637/jss.v100.i19."
  )
  model_description <- data.frame(model=names(models),
                                  inferred_parameter="theta (one scalar rate)",
                                  specification=c(rep("Fixed baseline; lambda and kappa calibrated to two quantiles", 3),
                                                  "Informative Gamma; shape and rate calibrated to two quantiles",
                                                  "Informative Lognormal; mu and sigma calibrated to two quantiles",
                                                  "Informative Weibull; shape and scale calibrated to two quantiles",
                                                  "Informative Gamma + weak Gamma; fixed weight 0.5",
                                                  "Same components; data-dependent SAM weight (empirical Bayes)"),
                                  comparison=c(rep("Same two quantiles", 6),
                                               "Robustification: marginal mixture quantiles are not matched",
                                               "Adaptation: data-dependent prior; marginal quantiles are not matched"))
  # Check grid results independently with adaptive quadrature.
  
  reference_integral <- function(m, a, b, x=cfg$threshold) {
    ob <- grid(m, a, b, N=cfg$grid_points)
    z0 <- min(ob$z)
    z1 <- max(ob$z)
    
    f <- function(z) exp(ob$lf(z)-ob$peak)
    
    ii <- function(fun, l=z0, h=z1) {
      if(h<=l) return(0)
      integrate(fun, l, h, rel.tol=1e-9, abs.tol=1e-11, subdivisions=2000)$value
    }
    den <- ii(f)
    
    F <- function(t) {
      zz <- log(t)
      if(zz<=z0) return(0)
      if(zz>=z1) return(1)
      ii(f, h=zz)/den
    }
    
    Q <- function(p) exp(uniroot(function(z)F(exp(z))-p, c(z0, z1), tol=1e-9)$root)
    c(mean=ii(function(z)exp(z)*f(z))/den,
      lower=Q(.025), upper=Q(.975), tail=1-F(x), F_x=F(x))
  }
  cat("1/6 Calibration and numerical checks\n")
  prior_table <- list()
  
  for(m in fixed_models) {
    cat("  Prior check:", m$name, "\n")
    flush.console()
    ob <- if(m$type=="mixture") NULL else grid(m, 0, 0)
    
    F <- function(t) {
      ob0 <- reference_integral(m, 0, 0, x=t)
      ob0["F_x"]
    }
    if(m$type=="mixture") {
      
      Q <- function(p) uniroot(function(t)
        cfg$mixture_weight*pgamma(t, gg[1], rate=gg[2])+
          (1-cfg$mixture_weight)*pgamma(t, a0, rate=b0)-p,
        c(1e-15, 100), tol=1e-11)$root
    } else Q <- m$quantile
    qs <- if(m$type=="mixture") vapply(c(.025, .5, .95, .975), Q, 0.0) else
      Q(c(.025, .5, .95, .975))
    prior_table[[m$name]] <- data.frame(model=m$name,
                                        lambda=if(is.null(m$lambda))NA_real_ else m$lambda,
                                        kappa=if(is.null(m$kappa))NA_real_ else m$kappa,
                                        q025=qs[1], q50=qs[2], q95=qs[3], q975=qs[4],
                                        prob_above_threshold=if(m$type=="mixture")
                                          cfg$mixture_weight*pgamma(cfg$threshold, gg[1], rate=gg[2], lower.tail=FALSE)+
                                          (1-cfg$mixture_weight)*pgamma(cfg$threshold, a0, rate=b0, lower.tail=FALSE)
                                        else 1-F(cfg$threshold))
  }
  prior_table <- bind_rows(prior_table)
  validations <- list()
  validation_models <- if (mode == "test") list(proposed_models[[1]], gam) else
    c(proposed_models, list(gam, ln, wei))
  validation_counts <- if (mode == "test") 0 else c(0, 2, 40)
  
  for(m in validation_models) for(a in validation_counts) {
    cat("  Numerical check:", m$name, "| count:", a, "\n")
    flush.console()
    b <- if(a==0)5 else a/.2
    num <- numerical(m, a, b)
    ref <- reference_integral(m, a, b)
    fields <- c("mean", "lower", "upper", "tail")
    errors <- abs(num[fields]-ref[fields])
    validations[[length(validations)+1]] <- data.frame(model=m$name, a=a, b=b,
                                                       mean_error=errors[1], q025_error=errors[2], q975_error=errors[3], tail_error=errors[4])
    if(max(errors)>2e-4) stop("Independent validation failed: ", m$name)
    if(m$type=="gamma") {
      exact <- posterior_gamma(m$a+a, m$b+b)
      if(max(abs(num[fields]-exact[fields]))>2e-4)
        stop("Quadrature failed to reproduce the analytical Gamma posterior")
    }
  }
  # Regression check for the reported zero-event censoring failure.
  for (m in proposed_models) {
    cat("  Zero-event regression check:", m$name, "\n")
    flush.console()
    a <- 0
    b <- 10.670861270364
    computed <- numerical(m, a, b)
    reference <- reference_integral(m, a, b)
    fields <- c("mean", "lower", "upper", "tail")
    errors <- abs(computed[fields] - reference[fields])
    
    if (max(errors) > 2e-4) {
      stop("Zero-event regression check failed: ", m$name)
    }
    
    validations[[length(validations) + 1L]] <- data.frame(
      model = m$name, a = a, b = b,
      mean_error = errors[1], q025_error = errors[2],
      q975_error = errors[3], tail_error = errors[4]
    )
  }
  
  validations <- bind_rows(validations)
  # RMSE MCSE uses the delta method; other means use sd/sqrt(M).
  
  summarize <- function(d, theta) {
    e <- d$mean-theta
    se <- e^2
    covered <- as.numeric(d$lower<=theta & theta<=d$upper)
    lengths <- d$upper-d$lower
    rm <- sqrt(mean(se))
    bs <- (d$tail-as.numeric(theta>cfg$threshold))^2
    signal <- as.numeric(d$tail>.95)
    lp <- -log(pmax(d$future_count_probability, .Machine$double.xmin))
    c(M=nrow(d), bias=mean(e), mcse_bias=mcse(e), RMSE=rm,
      mcse_RMSE=if(rm>0)mcse(se)/(2*rm) else 0,
      coverage=mean(covered), mcse_coverage=sqrt(mean(covered)*(1-mean(covered))/length(covered)),
      length=mean(lengths), mcse_length=mcse(lengths),
      mean_tail=mean(d$tail), mcse_tail=mcse(d$tail),
      threshold_Brier=mean(bs), mcse_Brier=mcse(bs),
      prob_signal=mean(signal), mcse_signal=sqrt(mean(signal)*(1-mean(signal))/length(signal)),
      predictive_NLL=mean(lp), mcse_NLL=mcse(lp),
      survival_error=mean(d$survival-exp(-cfg$horizon*theta)),
      mcse_surv=mcse(d$survival-exp(-cfg$horizon*theta)))
  }
  
  future_probability <- function(m, a, b, y) {
    ef <- cfg$future_exposure
    if(m$type=="gamma") return(dnbinom(y, size=m$a+a, prob=(m$b+b)/(m$b+b+ef)))
    if(m$type %in% c("mixture", "sam")) {
      p <- posterior_mix(m, a, b)["posterior_weight"]
      return(sum(c(p, 1-p)*dnbinom(y, size=m$a+a, prob=(m$b+b)/(m$b+b+ef))))
    }
    ob <- numerical(m, a, b, return_grid=TRUE)
    
    if (!is.null(ob$expectation)) {
      return(ob$expectation(function(theta) dpois(y, ef * theta)))
    }
    
    trap(ob$z, dpois(y, ef*exp(ob$z))*ob$density_z)
  }
  summaries <- list()
  paired_differences <- list()
  raw_results <- list()
  batch_cost <- list()
  # Cache Poisson posteriors by model, exposure and observed count.
  cache <- new.env(parent=emptyenv())
  cat("2/6 Poisson Monte Carlo (shared data across models)\n")
  cached <- stage_read("poisson")
  if (!is.null(cached)) {
    list2env(cached, envir=environment())
  } else {
  set.seed(cfg$seed+1L)
  
  for(E in cfg$exposures) for(theta in cfg$true_rates) {
    cat("Exposure =", E, "; true rate =", theta, "\n")
    y <- rpois(cfg$M, E*theta)
    yf <- rpois(cfg$M, cfg$future_exposure*theta)
    by_model <- list()
    
    for(m in models) {
      timer <- proc.time()[3]
      totals <- sort(unique(y))
      progress <- make_progress(paste(m$name, "posteriors"), length(totals), 25L)
      progress(0L)
      
      tab <- lapply(seq_along(totals), function(index) {
        a <- totals[index]
        progress(index)
        key <- paste(m$name, E, a, sep="|")
        if(!exists(key, cache, inherits=FALSE)) assign(key, posterior(m, a, E), cache)
        get(key, cache, inherits=FALSE)
      })
      tab <- do.call(rbind, tab)
      d <- as.data.frame(tab[match(y, totals), , drop=FALSE])
      ch <- paste(y, yf, sep=":")
      un <- !duplicated(ch)
      prediction_indices <- which(un)
      progress <- make_progress(paste(m$name, "predictions"), length(prediction_indices), 25L)
      progress(0L)
      
      pr <- vapply(seq_along(prediction_indices), function(index) {
        progress(index)
        i <- prediction_indices[index]
        future_probability(m, y[i], E, yf[i])
      }, 0.0)
      d$future_count_probability <- pr[match(ch, ch[un])]
      d$rep <- seq_len(cfg$M)
      d$y <- y
      d$y_futuro <- yf
      by_model[[m$name]] <- d
      summaries[[length(summaries)+1]] <- data.frame(exposure=E, theta=theta,
                                                     expected_events=E*theta, model=m$name, as.list(summarize(d, theta)), check.names=FALSE)
      batch_cost[[length(batch_cost)+1]] <- data.frame(exposure=E, theta=theta, model=m$name,
                                                       seconds=unname(proc.time()[3]-timer), distinct_totals=length(totals),
                                                       note="Cached batch time; not a per-posterior benchmark")
    }
    ref <- by_model[[gam$name]]
    
    for(nm in names(by_model)) {
      d <- by_model[[nm]]
      delta_mse <- (d$mean-theta)^2-(ref$mean-theta)^2
      delta_brier <- (d$tail-as.numeric(theta>cfg$threshold))^2-
        (ref$tail-as.numeric(theta>cfg$threshold))^2
      paired_differences[[length(paired_differences)+1]] <- data.frame(exposure=E, theta=theta, model=nm,
                                                                       delta_MSE=mean(delta_mse), mcse_delta_MSE=mcse(delta_mse),
                                                                       delta_Brier=mean(delta_brier), mcse_delta_Brier=mcse(delta_brier))
    }
    raw_results[[paste(E, theta, sep="|")]] <- by_model
  }
  summaries <- bind_rows(summaries)
  paired_differences <- bind_rows(paired_differences)
  batch_cost <- bind_rows(batch_cost)
    stage_save("poisson", c("summaries", "paired_differences", "raw_results", "batch_cost"), environment())
  }

  cat("3/6 SBC for each QROG baseline\n")
  cached <- stage_read("sbc")
  if (!is.null(cached)) {
    list2env(cached, envir=environment())
  } else {
  set.seed(cfg$seed+2L)
  # SBC uses numerical posterior inversion; retain the discrete rank probabilities.
  sbc <- list()
  sbc_tables <- list()
  sbc_validation <- list()
  
  for(m in proposed_models) {
    cat("SBC:", m$name, "\n")
    ranks <- integer(cfg$R_sbc)
    pit <- numeric(cfg$R_sbc)
    progress <- make_progress(paste(m$name, "SBC"), cfg$R_sbc)
    progress(0L)
    
    for(i in seq_len(cfg$R_sbc)) {
      progress(i)
      theta <- m$sample(1)
      a <- rpois(1, 20*theta)
      ob <- numerical(m, a, 20, return_grid=TRUE)
      draws <- ob$q(runif(cfg$L_sbc))
      ranks[i] <- sum(draws<theta)
      pit[i] <- ob$F(theta)
      if(i<=min(20, cfg$R_sbc)) {
        ref <- reference_integral(m, a, 20, x=theta)
        sbc_validation[[length(sbc_validation)+1]] <- data.frame(model=m$name,
                                                                 rep=i, grid_PIT=pit[i], PIT_integrate=ref["F_x"], error=abs(pit[i]-ref["F_x"]))
        if(abs(pit[i]-ref["F_x"])>2e-4) stop("SBC PIT check failed")
      }
    }
    K <- cfg$L_sbc+1L
    bins <- floor((0:cfg$L_sbc)*10/K)+1L
    observed <- tabulate(bins[ranks+1L], nbins=10)
    expected <- cfg$R_sbc*tabulate(bins, nbins=10)/K
    test_statistic <- sum((observed-expected)^2/expected)
    full_counts <- tabulate(ranks+1L, nbins=K)
    full_chi <- sum((full_counts-cfg$R_sbc/K)^2/(cfg$R_sbc/K))
    sbc[[m$name]] <- list(ranks=ranks, PIT=pit,
                          rank_mean=mean(ranks), theoretical_mean=cfg$L_sbc/2,
                          rank_variance=var(ranks), theoretical_variance=cfg$L_sbc*(cfg$L_sbc+2)/12,
                          chi2_10=test_statistic, p_10=if (mode == "test") NA_real_ else pchisq(test_statistic, 9, lower.tail=FALSE),
                          chi2_full_ranks=full_chi,
                          p_full_ranks=if (mode == "test") NA_real_ else pchisq(full_chi, K-1, lower.tail=FALSE),
                          KS_PIT=ks.test(pit, "punif", exact=FALSE))
    sbc_tables[[m$name]] <- data.frame(model=m$name, bin=1:10,
                                       observed=observed, expected=expected)
  }
  sbc_tables <- bind_rows(sbc_tables)
  sbc_validation <- bind_rows(sbc_validation)
  
    stage_save("sbc", c("sbc", "sbc_tables", "sbc_validation"), environment())
  }

  cat("4/6 Independent censoring: repeated experiments\n")
  cached <- stage_read("censoring")
  if (!is.null(cached)) {
    list2env(cached, envir=environment())
  } else {
  set.seed(cfg$seed+3L)
  censoring_summary <- list()
  censoring_raw <- list()
  
  for(n in cfg$censoring_sizes) for(frac in cfg$censoring_fractions)
    
    for(theta in cfg$censoring_rates) {
      cat("n =", n, "; expected censoring =", frac, "; rate =", theta, "\n")
      # Set independent censoring to the target fraction; omit its rate from inference.
      gamma_c <- theta*frac/(1-frac)
      data <- replicate(cfg$M_censoring, {
        tt <- rexp(n, theta)
        cc <- rexp(n, gamma_c)
        c(a=sum(tt<=cc), b=sum(pmin(tt, cc)))
      })
      
      for(m in models) {
        progress <- make_progress(paste(m$name, "censoring posteriors"), cfg$M_censoring)
        progress(0L)
        
        d <- bind_rows(lapply(seq_len(cfg$M_censoring), function(i) {
          progress(i)
          fitted <- posterior(m, data["a", i], data["b", i])
          as.data.frame(as.list(fitted))
        }))
        d$future_count_probability <- 1
        su <- summarize(d, theta)
        su <- su[!names(su)%in%c("predictive_NLL", "mcse_NLL")]
        censoring_summary[[length(censoring_summary)+1]] <- data.frame(n=n, censoring=frac, theta=theta,
                                                                       mean_events=mean(data["a", ]), mean_exposure=mean(data["b", ]),
                                                                       model=m$name, as.list(su), check.names=FALSE)
        censoring_raw[[paste(n, frac, theta, m$name, sep="|")]] <- list(statistics=data, posterior=d)
      }
    }
  censoring_summary <- bind_rows(censoring_summary)
  
    stage_save("censoring", c("censoring_summary", "censoring_raw"), environment())
  }

  cat("5/6 Special-case rejection efficiency and benchmarks\n")
  cached <- stage_read("rejection")
  if (!is.null(cached)) {
    list2env(cached, envir=environment())
  } else {
  set.seed(cfg$seed+4L)
  m0 <- proposed_models[[1]]
  s <- cfg$median
  # Exact rejection applies to baseline LL1 on the scaled parameter theta/s.
  
  rejection_sample <- function(a, b, N) {
    lam <- m0$lambda
    k <- m0$kappa
    rate <- b+lam*(1+k)/s
    
    logA <- function(t) log(qp(t/s, k))-log1p(k)-lam*k/(1+t/s)
    prob <- integrate(function(u)exp(logA(qgamma(u, a+1, rate=rate))), 0, 1,
                      rel.tol=1e-8, subdivisions=2000)$value
    samples <- numeric()
    n_proposals <- accepted <- 0L
    t0 <- proc.time()[3]
    while(length(samples)<N) {
      # Count all accepted proposals, including surplus draws in the last batch.
      batch <- min(200000L, max(1000L, ceiling((N-length(samples))/max(prob, 1e-5)*1.1)))
      t <- rgamma(batch, a+1, rate=rate)
      ok <- log(runif(batch))<=logA(t)
      n_proposals <- n_proposals+batch
      accepted <- accepted+sum(ok)
      samples <- c(samples, t[ok])
      if(n_proposals>5e7) stop("Rejection sampler exceeded its proposal limit")
    }
    elapsed <- proc.time()[3]-t0
    list(draws=samples[seq_len(N)], summary=data.frame(a=a, b=b,
                                                       theoretical_probability=prob, empirical_probability=accepted/n_proposals,
                                                       mcse_acceptance=sqrt(prob*(1-prob)/n_proposals), proposals=n_proposals,
                                                       total_accepted=accepted, retained=N, seconds=unname(elapsed)))
  }
  efficiency <- list()
  benchmarks <- list()
  
  for(E in cfg$exposures) for(theta in cfg$benchmark_rates) {
    a <- round(E*theta)
    r <- rejection_sample(a, E, cfg$N_rejection)
    ex <- posterior(m0, a, E)
    r$summary$mean_error <- mean(r$draws)-ex["mean"]
    r$summary$mcse_mean <- sd(r$draws)/sqrt(length(r$draws))
    r$summary$representative_rate <- theta
    efficiency[[length(efficiency)+1]] <- r$summary
    # Quadrature and sampling solve different tasks; report their costs separately.
    t1 <- proc.time()[3]
    ref <- reference_integral(m0, a, E)
    elapsed_quadrature <- proc.time()[3]-t1
    t1 <- proc.time()[3]
    ob <- numerical(m0, a, E, return_grid=TRUE)
    draws <- ob$q(runif(cfg$N_rejection))
    elapsed_inversion <- proc.time()[3]-t1
    benchmarks[[length(benchmarks)+1]] <- data.frame(a=a, b=E,
                                                     representative_rate=theta, adaptive_quadrature_seconds=unname(elapsed_quadrature),
                                                     numerical_inversion_with_setup_seconds=unname(elapsed_inversion),
                                                     rejection_N_seconds=r$summary$seconds, N=cfg$N_rejection,
                                                     grid_quantile_error=max(abs(ob$summary[c("lower", "upper")]-ref[c("lower", "upper")])),
                                                     inversion_mean_error=mean(draws)-ex["mean"], mcse_inversion_mean=sd(draws)/sqrt(length(draws)))
  }
  efficiency <- bind_rows(efficiency)
  benchmarks <- bind_rows(benchmarks)
  
    stage_save("rejection", c("efficiency", "benchmarks"), environment())
  }

  cat("6/6 Saving three output files\n")
  flush.console()
  notes <- c(
    paste("Run mode:", mode, "- test outputs are workflow checks only."),
    "Illustrative prior inputs; no clinical elicitation or historical meta-analysis is claimed.",
    "QROG, Gamma, Lognormal and Weibull share two quantiles; robust mixtures do not.",
    "SAM is an empirical-Bayes adaptation of the published likelihood-ratio rule.",
    "SBC uses numerical inversion; exact rejection applies only to the special baseline.",
    "MCSEs quantify simulation uncertainty; benchmark tasks and additional mixture settings differ."
  )
  results <- list(configuration=cfg, models=model_description, priors=prior_table,
                  validations=validations, simulation=summaries, paired_differences=paired_differences,
                  raw_MC=raw_results, batch_times=batch_cost, SBC=sbc, SBC_histograms=sbc_tables,
                  SBC_validation=sbc_validation, censoring=censoring_summary, raw_censoring=censoring_raw,
                  rejection_efficiency=efficiency, benchmarks=benchmarks,
                  grid_diagnostics=list(refinements=n_refinements, adaptive_fallbacks=n_adaptive_fallbacks, largest_initial_difference=max_difference),
                  notes=notes, references=references, warnings_log=warnings_log, sessionInfo=sessionInfo(),
                  total_seconds=unname(proc.time()[3]-start_time))
  export_validation(results, output_dir)
  invisible(results)
}

# Run one of these commands in the RStudio console after loading the definitions.
# run_qrog("test") or run_qrog("execution")

# Rebuild prior curves from saved calibration settings; no simulations are repeated.
validation_prior_curves <- function(results) {
  cfg <- results$configuration
  theta <- exp(seq(log(.0001), log(3), length.out = 2000))
  curves <- list()
  
  for (index in seq_len(3L)) {
    calibration <- results$priors[index, ]
    
    if (index <= 2L) {
      beta <- cfg$baseline_shapes[index]
      odds <- (theta / cfg$median)^beta
      log_jacobian <- log(beta / cfg$median) +
        (beta - 1) * log(theta / cfg$median)
    } else {
      z <- (log(theta) - log(cfg$median)) / cfg$sdlog_baseline
      log_survival <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
      odds <- exp(pnorm(z, log.p = TRUE) - log_survival)
      log_jacobian <- dnorm(z, log = TRUE) - log(theta) -
        log(cfg$sdlog_baseline) - 2 * log_survival
    }
    
    kappa <- calibration$kappa
    lambda <- calibration$lambda
    log_density <- log(lambda) + log(1 + kappa - kappa / (1 + odds)^2) +
      log_jacobian - lambda * (odds + kappa * odds^2 / (1 + odds))
    
    curves[[index]] <- data.frame(
      theta = theta, density = exp(log_density), model = calibration$model
    )
  }
  
  probabilities <- cfg$probabilities
  target_ratio <- cfg$upper / cfg$median
  log_shape <- uniroot(function(value) {
    shape <- exp(value)
    log(qgamma(probabilities[2], shape) / qgamma(probabilities[1], shape)) -
      log(target_ratio)
  }, c(log(.02), log(10000)), tol = 1e-11)$root
  
  shape <- exp(log_shape)
  rate <- qgamma(probabilities[1], shape) / cfg$median
  gamma_density <- dgamma(theta, shape, rate = rate)
  sigma <- log(target_ratio) / qnorm(probabilities[2])
  weibull_shape <- log(log1p(-probabilities[2]) / log1p(-probabilities[1])) /
    log(target_ratio)
  weibull_scale <- cfg$median / (-log1p(-probabilities[1]))^(1 / weibull_shape)
  weak_shape <- .25
  weak_rate <- weak_shape / (shape / rate)
  
  curves[[4L]] <- data.frame(theta, density = gamma_density, model = "Gamma-Poisson")
  curves[[5L]] <- data.frame(
    theta, density = dlnorm(theta, log(cfg$median), sigma), model = "Lognormal-Poisson"
  )
  curves[[6L]] <- data.frame(
    theta, density = dweibull(theta, weibull_shape, weibull_scale), model = "Weibull-Poisson"
  )
  curves[[7L]] <- data.frame(
    theta,
    density = cfg$mixture_weight * gamma_density +
      (1 - cfg$mixture_weight) * dgamma(theta, weak_shape, rate = weak_rate),
    model = "Robust-mixture"
  )
  
  do.call(rbind, curves)
}

# Model styles remain identical across all PDF pages, including grayscale line patterns.
plot_validation <- function(results, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('Install first: install.packages("ggplot2")')
  }
  
  cfg <- results$configuration
  model_order <- c(
    "QROG-LL1", "QROG-LL0.85", "QROG-LN2", "Gamma-Poisson",
    "Lognormal-Poisson", "Weibull-Poisson", "Robust-mixture", "SAM-Poisson"
  )
  colors <- setNames(c(
    "#D55E00", "#0072B2", "#009E73", "#000000",
    "#CC79A7", "#E69F00", "#56B4E9", "#666666"
  ), model_order)
  line_types <- setNames(c(
    "solid", "longdash", "dotdash", "dashed",
    "dotted", "twodash", "4313", "3313"
  ), model_order)
  shapes <- setNames(c(16, 17, 15, 1, 2, 0, 18, 8), model_order)
  
  figure_theme <- ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = .35),
      panel.border = ggplot2::element_rect(color = "grey35", linewidth = .6),
      axis.title = ggplot2::element_text(size = 15),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", linewidth = .5),
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.width = grid::unit(1.5, "cm"),
      legend.key.height = grid::unit(.55, "cm"),
      panel.spacing = grid::unit(.9, "lines"),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
  
  model_scales <- function(points = TRUE) {
    scales <- list(
      ggplot2::scale_color_manual(values = colors, breaks = model_order, name = NULL),
      ggplot2::scale_linetype_manual(values = line_types, breaks = model_order, name = NULL)
    )
    
    if (points) {
      scales <- c(scales, list(
        ggplot2::scale_shape_manual(values = shapes, breaks = model_order, name = NULL)
      ))
    }
    
    scales
  }
  
  legend_guides <- ggplot2::guides(
    color = ggplot2::guide_legend(ncol = 3, byrow = TRUE),
    linetype = ggplot2::guide_legend(ncol = 3, byrow = TRUE),
    shape = ggplot2::guide_legend(ncol = 3, byrow = TRUE)
  )
  
  prepare <- function(data) {
    data$model <- factor(data$model, levels = model_order)
    data[order(data$model), , drop = FALSE]
  }
  
  pdf_path <- file.path(output_dir, "simulation_figures.pdf")
  grDevices::pdf(pdf_path, width = 9.5, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  prior_data <- prepare(validation_prior_curves(results))
  print(
    ggplot2::ggplot(prior_data, ggplot2::aes(
      theta, density, color = model, linetype = model, group = model
    )) +
      ggplot2::geom_vline(xintercept = cfg$threshold, linetype = "dashed",
                          color = "grey45", linewidth = .7) +
      ggplot2::geom_line(linewidth = 1.15) +
      ggplot2::coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 10)) +
      model_scales(points = FALSE) + legend_guides + figure_theme +
      ggplot2::labs(title = "Prior distributions",
                    subtitle = "Common Q50 and Q95 for non-mixture priors; robust mixture not matched",
                    x = "Event rate", y = "Prior density")
  )
  
  metric_labels <- c(
    RMSE = "Root mean squared error", coverage = "95% interval coverage",
    mean_tail = "Posterior probability of exceeding the threshold",
    threshold_Brier = "Threshold quadratic risk", predictive_NLL = "Predictive negative log score"
  )
  error_columns <- c(
    RMSE = "mcse_RMSE", coverage = "mcse_coverage", mean_tail = "mcse_tail",
    threshold_Brier = "mcse_Brier", predictive_NLL = "mcse_NLL"
  )
  
  for (metric in names(metric_labels)) {
    data <- prepare(results$simulation)
    data$value <- data[[metric]]
    data$error <- data[[error_columns[[metric]]]]
    data$rate_panel <- paste("True rate =", data$theta)
    facet_scales <- if (metric %in% c("coverage", "mean_tail")) "fixed" else "free_y"
    
    figure <- ggplot2::ggplot(data, ggplot2::aes(
      exposure, value, color = model, linetype = model, shape = model, group = model
    )) +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = value - 1.96 * error, ymax = value + 1.96 * error
      ), width = 0, linewidth = .55, alpha = .65, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::geom_point(size = 2.6, stroke = .85) +
      ggplot2::scale_x_log10(breaks = cfg$exposures) +
      ggplot2::facet_wrap(~rate_panel, ncol = 2, scales = facet_scales) +
      model_scales() + legend_guides + figure_theme +
      ggplot2::labs(title = unname(metric_labels[[metric]]),
                    subtitle = "Error bars: +/- 1.96 Monte Carlo standard errors",
                    x = "Total exposure (log scale)", y = unname(metric_labels[[metric]]))
    
    if (metric == "coverage") {
      figure <- figure + ggplot2::geom_hline(
        yintercept = .95, linetype = "dashed", color = "grey35", linewidth = .65
      )
    }
    
    if (metric == "mean_tail") {
      figure <- figure + ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(y = "Pr(rate > 0.4 | data)")
    }
    
    print(figure)
  }
  
  paired <- prepare(results$paired_differences)
  paired$rate_panel <- paste("True rate =", paired$theta)
  print(
    ggplot2::ggplot(paired, ggplot2::aes(
      exposure, delta_MSE, color = model, linetype = model, shape = model, group = model
    )) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = .65) +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = delta_MSE - 1.96 * mcse_delta_MSE,
        ymax = delta_MSE + 1.96 * mcse_delta_MSE
      ), width = 0, linewidth = .55, alpha = .65, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 2.6, stroke = .85) +
      ggplot2::scale_x_log10(breaks = cfg$exposures) +
      ggplot2::facet_wrap(~rate_panel, ncol = 2, scales = "free_y") +
      model_scales() + legend_guides + figure_theme +
      ggplot2::labs(title = "Paired MSE differences relative to Gamma",
                    subtitle = "Negative values favor the comparator; bars: +/- 1.96 MCSE",
                    x = "Total exposure (log scale)", y = "MSE difference")
  )
  
  print(
    ggplot2::ggplot(results$SBC_histograms, ggplot2::aes(bin, observed)) +
      ggplot2::geom_col(fill = "#0072B2", color = "white", linewidth = .35) +
      ggplot2::geom_point(ggplot2::aes(y = expected), color = "#D55E00", size = 2.8) +
      ggplot2::facet_wrap(~model, nrow = 1) + figure_theme +
      ggplot2::scale_x_continuous(breaks = c(1, 5, 10)) +
      ggplot2::labs(title = "Simulation-based calibration",
                    subtitle = "Points indicate expected counts for discrete ranks",
                    x = "Rank group", y = "Replication count")
  )
  
  censored <- prepare(results$censoring)
  censored$panel <- paste0("Rate = ", censored$theta,
                           "; censoring = ", 100 * censored$censoring, "%")
  print(
    ggplot2::ggplot(censored, ggplot2::aes(
      n, coverage, color = model, linetype = model, shape = model, group = model
    )) +
      ggplot2::geom_hline(yintercept = .95, linetype = "dashed", linewidth = .65) +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = coverage - 1.96 * mcse_coverage, ymax = coverage + 1.96 * mcse_coverage
      ), width = 0, linewidth = .55, alpha = .65, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 2.6, stroke = .85) +
      ggplot2::scale_x_log10(breaks = cfg$censoring_sizes) +
      ggplot2::facet_wrap(~panel, ncol = 2) +
      model_scales() + legend_guides + figure_theme +
      ggplot2::labs(title = "Coverage under independent right censoring",
                    subtitle = "Error bars: +/- 1.96 Monte Carlo standard errors",
                    x = "Sample size (log scale)", y = "95% interval coverage")
  )
  
  print(
    ggplot2::ggplot(results$rejection_efficiency, ggplot2::aes(
      theoretical_probability, empirical_probability
    )) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                           linewidth = 1, color = "grey40") +
      ggplot2::geom_point(size = 3.2, shape = 21, fill = "#0072B2", color = "black") +
      ggplot2::coord_equal() + figure_theme +
      ggplot2::labs(title = "Exact rejection sampling",
                    x = "Theoretical acceptance probability", y = "Observed acceptance rate")
  )
  
  invisible(pdf_path)
}

# Export an existing RDS or a results list without rerunning any simulation.
export_validation <- function(results, output_dir = NULL) {
  if (is.character(results) && length(results) == 1L) {
    input_path <- results
    results <- readRDS(input_path)
    if (is.null(output_dir)) output_dir <- dirname(input_path)
  }
  
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "results", "simulation")
    if (identical(results$configuration$mode, "test")) {
      output_dir <- file.path(output_dir, "test")
    }
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(output_dir, "simulation_results.rds"), compress = "xz")
  
  report <- capture.output({
    cat("GENERAL QROG VALIDATION\n\n")
    print(results$configuration)
    cat("\nNOTES\n")
    cat(results$notes, sep = "\n")
    
    for (name in c(
      "models", "priors", "validations", "simulation", "paired_differences",
      "SBC_histograms", "SBC_validation", "censoring", "rejection_efficiency",
      "benchmarks", "grid_diagnostics"
    )) {
      cat("\n\n", toupper(name), "\n", sep = "")
      print(results[[name]], row.names = FALSE, digits = 7)
    }
    
    cat("\nSBC SUMMARY\n")
    for (name in names(results$SBC)) {
      cat(name, "\n")
      print(results$SBC[[name]][setdiff(names(results$SBC[[name]]), c("ranks", "PIT"))])
    }
    
    cat("\nREFERENCES\n")
    cat(results$references, sep = "\n")
    cat("\n\nSIMULATION TIME (s):", results$total_seconds, "\n")
    print(results$sessionInfo)
  })
  
  writeLines(report, file.path(output_dir, "simulation_report.txt"), useBytes = TRUE)
  plot_validation(results, output_dir)
  cat("\nFiles saved in:", normalizePath(output_dir), "\n")
  cat("simulation_figures.pdf | simulation_report.txt | simulation_results.rds\n")
  invisible(results)
}

