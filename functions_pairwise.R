# ------------------------------------------------------------------------------
# functions_pairwise.R
#
# Pairwise U-statistic utilities aligned with the corrected attached AoS article.
#
# Implements:
# - SPDV_n: sample pairwise difference variance
# - SPFD_n: sample pairwise fourth-order dispersion
# - Pairwise estimator of fourth central moment:
#   mu4_pairwise = 6*SPFD_n - 3*SPDV_n^2
# - Pairwise plug-in estimator of asymptotic variance of sqrt(n)(SPDV_n - sigma^2):
#   avar_spdv_pairwise = 6*SPFD_n - 4*SPDV_n^2
#
# Corrected identities used:
# SPDV_n = choose(n,2)^(-1) sum_{i<j} (X_i-X_j)^2 / 2 = sample variance
# SPFD_n = choose(n,2)^(-1) sum_{i<j} (X_i-X_j)^4 / 12
# mu4 = 6*theta4 - 3*sigma^4
# asymptotic var for sqrt(n)(SPDV_n - sigma^2) = mu4 - sigma^4 = 6*theta4 - 4*sigma^4
# psi2(x) = ((x-mu)^2 - sigma^2)/2
# tau4^2 = 4 Var(psi4(X))
# ------------------------------------------------------------------------------

# =========================
# Basic validation helpers
# =========================

.check_numeric_vector <- function(x, min_n = 1L, finite_only = TRUE) {
  if (!is.numeric(x) || is.null(x)) {
    stop("x must be a numeric vector.")
  }
  if (length(x) < min_n) {
    stop(sprintf("x must have length at least %d.", min_n))
  }
  if (finite_only && any(!is.finite(x))) {
    stop("x must contain only finite values.")
  }
  invisible(TRUE)
}

.check_scalar_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x != as.integer(x)) {
    stop(sprintf("%s must be a single positive integer.", name))
  }
  invisible(TRUE)
}

.check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number in (0,1).")
  }
  invisible(TRUE)
}

.check_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% choices)) {
    stop(sprintf("%s must be one of: %s", name, paste(choices, collapse = ", ")))
  }
  invisible(TRUE)
}

.binom2 <- function(n) {
  n * (n - 1) / 2
}

.z_alpha <- function(level = 0.95) {
  stats::qnorm((1 + level) / 2)
}

.empty_ci <- function(
  estimate = NA_real_,
  se = NA_real_,
  lower = NA_real_,
  upper = NA_real_,
  level = NA_real_,
  method = NA_character_,
  target = NA_character_,
  success = FALSE,
  failure_reason = "unknown",
  extras = list()
) {
  out <- list(
    estimate = estimate,
    se = se,
    lower = lower,
    upper = upper,
    level = level,
    method = method,
    target = target,
    success = success,
    failure_reason = failure_reason
  )
  if (length(extras) > 0L) out <- c(out, extras)
  out
}

is_valid_ci <- function(ci) {
  if (is.null(ci)) return(FALSE)
  required <- c("lower", "upper", "success")
  if (!all(required %in% names(ci))) return(FALSE)
  if (!isTRUE(ci$success)) return(FALSE)
  if (!is.finite(ci$lower) || !is.finite(ci$upper)) return(FALSE)
  if (ci$lower > ci$upper) return(FALSE)
  TRUE
}

safe_quantile <- function(x, probs = c(0.025, 0.975), na_rm = TRUE, type = 7) {
  if (is.null(x) || length(x) == 0L) {
    return(list(
      success = FALSE,
      value = rep(NA_real_, length(probs)),
      failure_reason = "empty_input"
    ))
  }

  if (na_rm) x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(list(
      success = FALSE,
      value = rep(NA_real_, length(probs)),
      failure_reason = "no_finite_values"
    ))
  }

  qv <- tryCatch(
    as.numeric(stats::quantile(x, probs = probs, na.rm = FALSE, type = type, names = FALSE)),
    error = function(e) rep(NA_real_, length(probs))
  )

  if (any(!is.finite(qv))) {
    return(list(
      success = FALSE,
      value = rep(NA_real_, length(probs)),
      failure_reason = "quantile_failure"
    ))
  }

  list(success = TRUE, value = qv, failure_reason = "ok")
}

boot_resample <- function(x) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  sample(x, size = length(x), replace = TRUE)
}

# =========================
# Pairwise core statistics
# =========================

spdv <- function(x) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  stats::var(x)
}

spfd <- function(x) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)

  n <- length(x)
  if (n < 2L) stop("x must have length at least 2.")

  s1 <- sum(x)
  s2 <- sum(x^2)
  s3 <- sum(x^3)
  s4 <- sum(x^4)

  sum_diff4 <- n * s4 - 4 * s1 * s3 + 3 * s2^2
  sum_pair_kernel <- sum_diff4 / 12

  sum_pair_kernel / .binom2(n)
}

mu4_pairwise <- function(x) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  6 * spfd(x) - 3 * spdv(x)^2
}

avar_spdv_pairwise <- function(x, truncate_nonnegative = TRUE) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  val <- 6 * spfd(x) - 4 * spdv(x)^2
  if (truncate_nonnegative) val <- max(val, 0)
  val
}

# ====================================
# Classical moment estimators / helpers
# ====================================

mu4_classical <- function(x) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  xc <- x - mean(x)
  mean(xc^4)
}

kurtosis_raw_classical <- function(x) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  mu4_classical(x) / stats::var(x)^2
}

kurtosis_excess_classical <- function(x) {
  kurtosis_raw_classical(x) - 3
}

theta4_from_mu4_sigma2 <- function(mu4, sigma2) {
  mu4 / 6 + sigma2^2 / 2
}

theta4_normal_var1 <- function() {
  theta4_from_mu4_sigma2(mu4 = 3, sigma2 = 1)
}

mu4_normal_var1 <- function() {
  3
}

theta4_std_t <- function(df) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
    stop("df must be a single finite number > 2.")
  }
  if (df <= 4) return(Inf)
  theta4_from_mu4_sigma2(mu4 = mu4_std_t(df), sigma2 = 1)
}

mu4_std_t <- function(df) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
    stop("df must be a single finite number > 2.")
  }
  if (df <= 4) return(Inf)
  3 * (df - 2) / (df - 4)
}

theta4_contam_norm <- function(eps = 0.10, sd_clean = 1, sd_contam = 10, standardize = FALSE) {
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0 || eps > 1) {
    stop("eps must be in [0,1].")
  }
  p1 <- 1 - eps
  p2 <- eps
  sigma2_raw <- p1 * sd_clean^2 + p2 * sd_contam^2
  mu4_raw <- 3 * (p1 * sd_clean^4 + p2 * sd_contam^4)

  if (standardize) {
    mu4_std <- mu4_raw / sigma2_raw^2
    return(theta4_from_mu4_sigma2(mu4 = mu4_std, sigma2 = 1))
  }

  theta4_from_mu4_sigma2(mu4 = mu4_raw, sigma2 = sigma2_raw)
}

mu4_contam_norm <- function(eps = 0.10, sd_clean = 1, sd_contam = 10, standardize = FALSE) {
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0 || eps > 1) {
    stop("eps must be in [0,1].")
  }
  p1 <- 1 - eps
  p2 <- eps
  sigma2_raw <- p1 * sd_clean^2 + p2 * sd_contam^2
  mu4_raw <- 3 * (p1 * sd_clean^4 + p2 * sd_contam^4)

  if (standardize) {
    return(mu4_raw / sigma2_raw^2)
  }

  mu4_raw
}

sigma2_contam_norm <- function(eps = 0.10, sd_clean = 1, sd_contam = 10, standardize = FALSE) {
  if (standardize) return(1)
  (1 - eps) * sd_clean^2 + eps * sd_contam^2
}

# ===========================
# Projection-based estimators
# ===========================

projection_spdv <- function(x) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  ((x - mean(x))^2 - spdv(x)) / 2
}

se_spdv_projection <- function(x, truncate_nonnegative = TRUE) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)

  psi_hat <- projection_spdv(x)
  val <- 4 * mean(psi_hat^2)

  if (truncate_nonnegative) val <- max(val, 0)
  sqrt(val / length(x))
}

se_spdv_asymptotic <- function(x, truncate_nonnegative = TRUE) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  avar_hat <- avar_spdv_pairwise(x, truncate_nonnegative = truncate_nonnegative)
  sqrt(avar_hat / length(x))
}

projection_spfd_plugin <- function(x) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)

  mu_hat <- mean(x)
  xc <- x - mu_hat
  mu2_hat <- mean(xc^2)
  mu3_hat <- mean(xc^3)
  mu4_hat <- mean(xc^4)

  xc^4 / 12 + (xc^2 * mu2_hat) / 2 - (xc * mu3_hat) / 3 - mu4_hat / 12 - mu2_hat^2 / 2
}

tau4_sq_plugin <- function(x, assume_symmetric = FALSE, truncate_nonnegative = TRUE) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)

  xc <- x - mean(x)
  mu2_hat <- mean(xc^2)
  mu3_hat <- mean(xc^3)
  mu4_hat <- mean(xc^4)
  mu5_hat <- mean(xc^5)
  mu6_hat <- mean(xc^6)
  mu8_hat <- mean(xc^8)

  if (assume_symmetric) {
    val <- mu8_hat / 36 +
      (2 * mu2_hat * mu6_hat) / 3 +
      mu2_hat^2 * mu4_hat -
      mu4_hat^2 / 9 -
      mu2_hat^4
  } else {
    val <- mu8_hat / 36 +
      (2 * mu2_hat * mu6_hat) / 3 -
      (2 * mu3_hat * mu5_hat) / 9 +
      mu2_hat^2 * mu4_hat -
      mu4_hat^2 / 9 -
      mu2_hat^4 -
      (4 * mu2_hat * mu3_hat^2) / 3
  }

  if (truncate_nonnegative) val <- max(val, 0)
  val
}

se_spfd_projection <- function(
  x,
  assume_symmetric = FALSE,
  method = c("master_formula", "empirical_projection"),
  truncate_nonnegative = TRUE
) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  method <- match.arg(method)

  tau4_sq_hat <- switch(
    method,
    master_formula = tau4_sq_plugin(
      x,
      assume_symmetric = assume_symmetric,
      truncate_nonnegative = truncate_nonnegative
    ),
    empirical_projection = {
      psi_hat <- projection_spfd_plugin(x)
      val <- 4 * mean(psi_hat^2)
      if (truncate_nonnegative) val <- max(val, 0)
      val
    }
  )

  sqrt(tau4_sq_hat / length(x))
}

# ====================
# Data-generating laws
# ====================

r_std_t <- function(n, df) {
  .check_scalar_positive_integer(n, "n")
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
    stop("df must be a single finite number > 2.")
  }
  stats::rt(n, df = df) * sqrt((df - 2) / df)
}

r_norm1 <- function(n) {
  .check_scalar_positive_integer(n, "n")
  stats::rnorm(n, mean = 0, sd = 1)
}

r_contam_norm <- function(n, eps = 0.10, sd_clean = 1, sd_contam = 10, standardize = FALSE) {
  .check_scalar_positive_integer(n, "n")
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0 || eps > 1) {
    stop("eps must be in [0,1].")
  }

  z <- stats::runif(n) < eps
  x <- stats::rnorm(n, mean = 0, sd = sd_clean)
  x[z] <- stats::rnorm(sum(z), mean = 0, sd = sd_contam)

  if (standardize) {
    sd_mix <- sqrt((1 - eps) * sd_clean^2 + eps * sd_contam^2)
    x <- x / sd_mix
  }

  x
}

# ==================================
# Distribution parsing / truth helpers
# ==================================

parse_dist_name <- function(dist_name) {
  if (!is.character(dist_name) || length(dist_name) != 1L || is.na(dist_name)) {
    stop("dist_name must be a single non-missing character string.")
  }

  if (dist_name == "normal") {
    return(list(family = "normal", label = "normal", df = NA_real_))
  }

  if (dist_name %in% c("contam", "contam_norm")) {
    return(list(family = "contam_norm", label = "contam", df = NA_real_))
  }

  if (grepl("^t[0-9]+(\\.[0-9]+)?$", dist_name)) {
    df <- as.numeric(sub("^t", "", dist_name))
    return(list(family = "t", label = dist_name, df = df))
  }

  stop("Unsupported dist_name. Use 'normal', 'contam', 'contam_norm', or labels like 't8', 't12.5'.")
}

truth_sigma2 <- function(
  dist_name = "normal",
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE
) {
  info <- parse_dist_name(dist_name)

  if (info$family == "normal") return(1)
  if (info$family == "t") return(1)
  if (info$family == "contam_norm") {
    return(sigma2_contam_norm(
      eps = eps, sd_clean = sd_clean, sd_contam = sd_contam,
      standardize = standardize_contam
    ))
  }

  stop("Unsupported family.")
}

truth_theta4 <- function(
  dist_name = "normal",
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE
) {
  info <- parse_dist_name(dist_name)

  if (info$family == "normal") return(theta4_normal_var1())
  if (info$family == "t") return(theta4_std_t(info$df))
  if (info$family == "contam_norm") {
    return(theta4_contam_norm(
      eps = eps, sd_clean = sd_clean, sd_contam = sd_contam,
      standardize = standardize_contam
    ))
  }

  stop("Unsupported family.")
}

truth_mu4 <- function(
  dist_name = "normal",
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE
) {
  info <- parse_dist_name(dist_name)

  if (info$family == "normal") return(mu4_normal_var1())
  if (info$family == "t") return(mu4_std_t(info$df))
  if (info$family == "contam_norm") {
    return(mu4_contam_norm(
      eps = eps, sd_clean = sd_clean, sd_contam = sd_contam,
      standardize = standardize_contam
    ))
  }

  stop("Unsupported family.")
}

simulate_one_sample <- function(
  n,
  dist_name = "normal",
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE
) {
  info <- parse_dist_name(dist_name)

  if (info$family == "normal") {
    return(r_norm1(n))
  }

  if (info$family == "t") {
    return(r_std_t(n, df = info$df))
  }

  if (info$family == "contam_norm") {
    return(r_contam_norm(
      n = n, eps = eps, sd_clean = sd_clean, sd_contam = sd_contam,
      standardize = standardize_contam
    ))
  }

  stop("Unsupported family.")
}

# ===========================
# Generic bootstrap machinery
# ===========================

.bootstrap_stat <- function(x, stat_fun, B = 1999L) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")

  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    xb <- boot_resample(x)
    out[b] <- tryCatch(stat_fun(xb), error = function(e) NA_real_)
  }
  out
}

.bootstrap_studentized <- function(x, stat_fun, se_fun, B = 1999L) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")

  tstar <- rep(NA_real_, B)
  stat_star <- rep(NA_real_, B)
  se_star <- rep(NA_real_, B)

  stat_hat <- tryCatch(stat_fun(x), error = function(e) NA_real_)
  se_hat <- tryCatch(se_fun(x), error = function(e) NA_real_)

  for (b in seq_len(B)) {
    xb <- boot_resample(x)
    sb <- tryCatch(stat_fun(xb), error = function(e) NA_real_)
    seb <- tryCatch(se_fun(xb), error = function(e) NA_real_)

    stat_star[b] <- sb
    se_star[b] <- seb

    if (is.finite(sb) && is.finite(seb) && seb > 0 && is.finite(stat_hat)) {
      tstar[b] <- (sb - stat_hat) / seb
    } else {
      tstar[b] <- NA_real_
    }
  }

  list(
    stat_hat = stat_hat,
    se_hat = se_hat,
    stat_star = stat_star,
    se_star = se_star,
    t_star = tstar
  )
}

# ===========================================
# Confidence intervals for sigma^2 via SPDV_n
# ===========================================

ci_sigma2_normal <- function(
  x,
  level = 0.95,
  se_method = c("pairwise_plugin", "empirical_projection")
) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  .check_level(level)
  se_method <- match.arg(se_method)

  est <- tryCatch(spdv(x), error = function(e) NA_real_)
  se <- tryCatch(
    switch(
      se_method,
      pairwise_plugin = se_spdv_asymptotic(x),
      empirical_projection = se_spdv_projection(x)
    ),
    error = function(e) NA_real_
  )

  if (!is.finite(est)) {
    return(.empty_ci(
      estimate = est, se = se, level = level, method = "normal",
      target = "sigma2", failure_reason = "nonfinite_estimate",
      extras = list(se_method = se_method)
    ))
  }

  if (!is.finite(se) || se <= 0) {
    return(.empty_ci(
      estimate = est, se = se, level = level, method = "normal",
      target = "sigma2", failure_reason = "nonpositive_se",
      extras = list(se_method = se_method)
    ))
  }

  z <- .z_alpha(level)
  .empty_ci(
    estimate = est,
    se = se,
    lower = est - z * se,
    upper = est + z * se,
    level = level,
    method = "normal",
    target = "sigma2",
    success = TRUE,
    failure_reason = "ok",
    extras = list(se_method = se_method)
  )
}

ci_sigma2_percentile <- function(x, B = 1999L, level = 0.95) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")
  .check_level(level)

  est <- tryCatch(spdv(x), error = function(e) NA_real_)
  if (!is.finite(est)) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "sigma2", failure_reason = "nonfinite_estimate"
    ))
  }

  boot_stats <- .bootstrap_stat(x, spdv, B = B)
  boot_finite <- boot_stats[is.finite(boot_stats)]

  if (length(boot_finite) < max(20L, ceiling(0.5 * B))) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "sigma2", failure_reason = "too_few_finite_bootstrap_reps",
      extras = list(bootstrap_reps_finite = length(boot_finite))
    ))
  }

  alpha <- 1 - level
  qout <- safe_quantile(boot_finite, probs = c(alpha / 2, 1 - alpha / 2))
  if (!qout$success) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "sigma2", failure_reason = "bootstrap_quantile_failure",
      extras = list(bootstrap_reps_finite = length(boot_finite))
    ))
  }

  .empty_ci(
    estimate = est,
    se = stats::sd(boot_finite),
    lower = qout$value[1L],
    upper = qout$value[2L],
    level = level,
    method = "percentile",
    target = "sigma2",
    success = TRUE,
    failure_reason = "ok",
    extras = list(bootstrap_reps_finite = length(boot_finite))
  )
}

ci_sigma2_studentized <- function(
  x,
  B = 1999L,
  level = 0.95,
  se_method = c("pairwise_plugin", "empirical_projection")
) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")
  .check_level(level)
  se_method <- match.arg(se_method)

  stat_fun <- spdv
  se_fun <- switch(
    se_method,
    pairwise_plugin = se_spdv_asymptotic,
    empirical_projection = se_spdv_projection
  )

  bs <- .bootstrap_studentized(x, stat_fun = stat_fun, se_fun = se_fun, B = B)

  if (!is.finite(bs$stat_hat)) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "sigma2",
      failure_reason = "nonfinite_estimate",
      extras = list(se_method = se_method)
    ))
  }

  if (!is.finite(bs$se_hat) || bs$se_hat <= 0) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "sigma2",
      failure_reason = "nonpositive_se",
      extras = list(se_method = se_method)
    ))
  }

  t_finite <- bs$t_star[is.finite(bs$t_star)]
  if (length(t_finite) < max(20L, ceiling(0.5 * B))) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "sigma2",
      failure_reason = "too_few_finite_bootstrap_reps",
      extras = list(
        se_method = se_method,
        bootstrap_reps_finite = length(t_finite),
        bootstrap_se_finite = sum(is.finite(bs$se_star))
      )
    ))
  }

  alpha <- 1 - level
  qout <- safe_quantile(t_finite, probs = c(1 - alpha / 2, alpha / 2))
  if (!qout$success) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "sigma2",
      failure_reason = "bootstrap_quantile_failure",
      extras = list(
        se_method = se_method,
        bootstrap_reps_finite = length(t_finite)
      )
    ))
  }

  lower <- bs$stat_hat - qout$value[1L] * bs$se_hat
  upper <- bs$stat_hat - qout$value[2L] * bs$se_hat

  .empty_ci(
    estimate = bs$stat_hat,
    se = bs$se_hat,
    lower = lower,
    upper = upper,
    level = level,
    method = "studentized",
    target = "sigma2",
    success = TRUE,
    failure_reason = "ok",
    extras = list(
      se_method = se_method,
      bootstrap_reps_finite = length(t_finite),
      bootstrap_sd_t = stats::sd(t_finite)
    )
  )
}

# =====================================================
# Confidence intervals for theta4 via SPFD_n (quartic)
# =====================================================

ci_theta4_normal <- function(
  x,
  level = 0.95,
  se_method = c("master_formula", "empirical_projection"),
  assume_symmetric = FALSE
) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  .check_level(level)
  se_method <- match.arg(se_method)

  est <- tryCatch(spfd(x), error = function(e) NA_real_)
  se <- tryCatch(
    se_spfd_projection(
      x,
      assume_symmetric = assume_symmetric,
      method = se_method
    ),
    error = function(e) NA_real_
  )

  if (!is.finite(est)) {
    return(.empty_ci(
      estimate = est, se = se, level = level, method = "normal",
      target = "theta4", failure_reason = "nonfinite_estimate",
      extras = list(se_method = se_method, assume_symmetric = assume_symmetric)
    ))
  }

  if (!is.finite(se) || se <= 0) {
    return(.empty_ci(
      estimate = est, se = se, level = level, method = "normal",
      target = "theta4", failure_reason = "nonpositive_se",
      extras = list(se_method = se_method, assume_symmetric = assume_symmetric)
    ))
  }

  z <- .z_alpha(level)
  .empty_ci(
    estimate = est,
    se = se,
    lower = est - z * se,
    upper = est + z * se,
    level = level,
    method = "normal",
    target = "theta4",
    success = TRUE,
    failure_reason = "ok",
    extras = list(se_method = se_method, assume_symmetric = assume_symmetric)
  )
}

ci_theta4_percentile <- function(x, B = 1999L, level = 0.95) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")
  .check_level(level)

  est <- tryCatch(spfd(x), error = function(e) NA_real_)
  if (!is.finite(est)) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "theta4", failure_reason = "nonfinite_estimate"
    ))
  }

  boot_stats <- .bootstrap_stat(x, spfd, B = B)
  boot_finite <- boot_stats[is.finite(boot_stats)]

  if (length(boot_finite) < max(20L, ceiling(0.5 * B))) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "theta4", failure_reason = "too_few_finite_bootstrap_reps",
      extras = list(bootstrap_reps_finite = length(boot_finite))
    ))
  }

  alpha <- 1 - level
  qout <- safe_quantile(boot_finite, probs = c(alpha / 2, 1 - alpha / 2))
  if (!qout$success) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "theta4", failure_reason = "bootstrap_quantile_failure",
      extras = list(bootstrap_reps_finite = length(boot_finite))
    ))
  }

  .empty_ci(
    estimate = est,
    se = stats::sd(boot_finite),
    lower = qout$value[1L],
    upper = qout$value[2L],
    level = level,
    method = "percentile",
    target = "theta4",
    success = TRUE,
    failure_reason = "ok",
    extras = list(bootstrap_reps_finite = length(boot_finite))
  )
}

ci_theta4_studentized <- function(
  x,
  B = 1999L,
  level = 0.95,
  se_method = c("master_formula", "empirical_projection"),
  assume_symmetric = FALSE
) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")
  .check_level(level)
  se_method <- match.arg(se_method)

  se_fun <- function(z) {
    se_spfd_projection(
      z,
      assume_symmetric = assume_symmetric,
      method = se_method
    )
  }

  bs <- .bootstrap_studentized(x, stat_fun = spfd, se_fun = se_fun, B = B)

  if (!is.finite(bs$stat_hat)) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "theta4",
      failure_reason = "nonfinite_estimate",
      extras = list(se_method = se_method, assume_symmetric = assume_symmetric)
    ))
  }

  if (!is.finite(bs$se_hat) || bs$se_hat <= 0) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "theta4",
      failure_reason = "nonpositive_se",
      extras = list(se_method = se_method, assume_symmetric = assume_symmetric)
    ))
  }

  t_finite <- bs$t_star[is.finite(bs$t_star)]
  if (length(t_finite) < max(20L, ceiling(0.5 * B))) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "theta4",
      failure_reason = "too_few_finite_bootstrap_reps",
      extras = list(
        se_method = se_method,
        assume_symmetric = assume_symmetric,
        bootstrap_reps_finite = length(t_finite),
        bootstrap_se_finite = sum(is.finite(bs$se_star))
      )
    ))
  }

  alpha <- 1 - level
  qout <- safe_quantile(t_finite, probs = c(1 - alpha / 2, alpha / 2))
  if (!qout$success) {
    return(.empty_ci(
      estimate = bs$stat_hat, se = bs$se_hat, level = level,
      method = "studentized", target = "theta4",
      failure_reason = "bootstrap_quantile_failure",
      extras = list(
        se_method = se_method,
        assume_symmetric = assume_symmetric,
        bootstrap_reps_finite = length(t_finite)
      )
    ))
  }

  lower <- bs$stat_hat - qout$value[1L] * bs$se_hat
  upper <- bs$stat_hat - qout$value[2L] * bs$se_hat

  .empty_ci(
    estimate = bs$stat_hat,
    se = bs$se_hat,
    lower = lower,
    upper = upper,
    level = level,
    method = "studentized",
    target = "theta4",
    success = TRUE,
    failure_reason = "ok",
    extras = list(
      se_method = se_method,
      assume_symmetric = assume_symmetric,
      bootstrap_reps_finite = length(t_finite),
      bootstrap_sd_t = stats::sd(t_finite)
    )
  )
}

# ==========================================
# Fourth central moment intervals / comparison
# ==========================================

ci_mu4_pairwise_percentile <- function(x, B = 1999L, level = 0.95) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  .check_scalar_positive_integer(B, "B")
  .check_level(level)

  est <- tryCatch(mu4_pairwise(x), error = function(e) NA_real_)
  if (!is.finite(est)) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "mu4", failure_reason = "nonfinite_estimate"
    ))
  }

  boot_stats <- .bootstrap_stat(x, mu4_pairwise, B = B)
  boot_finite <- boot_stats[is.finite(boot_stats)]

  if (length(boot_finite) < max(20L, ceiling(0.5 * B))) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "mu4", failure_reason = "too_few_finite_bootstrap_reps"
    ))
  }

  alpha <- 1 - level
  qout <- safe_quantile(boot_finite, probs = c(alpha / 2, 1 - alpha / 2))
  if (!qout$success) {
    return(.empty_ci(
      estimate = est, level = level, method = "percentile",
      target = "mu4", failure_reason = "bootstrap_quantile_failure"
    ))
  }

  .empty_ci(
    estimate = est,
    se = stats::sd(boot_finite),
    lower = qout$value[1L],
    upper = qout$value[2L],
    level = level,
    method = "percentile",
    target = "mu4",
    success = TRUE,
    failure_reason = "ok"
  )
}

# =====================
# Monte Carlo wrappers
# =====================

compute_sigma2_intervals <- function(
  x,
  level = 0.95,
  B = 1999L,
  include_normal = TRUE,
  include_percentile = TRUE,
  include_studentized = TRUE,
  se_method = c("pairwise_plugin", "empirical_projection")
) {
  se_method <- match.arg(se_method)

  out <- list()
  if (include_normal) {
    out$normal <- ci_sigma2_normal(x, level = level, se_method = se_method)
  }
  if (include_percentile) {
    out$percentile <- ci_sigma2_percentile(x, B = B, level = level)
  }
  if (include_studentized) {
    out$studentized <- ci_sigma2_studentized(
      x, B = B, level = level, se_method = se_method
    )
  }
  out
}

compute_theta4_intervals <- function(
  x,
  level = 0.95,
  B = 1999L,
  include_normal = TRUE,
  include_percentile = TRUE,
  include_studentized = TRUE,
  se_method = c("master_formula", "empirical_projection"),
  assume_symmetric = FALSE
) {
  se_method <- match.arg(se_method)

  out <- list()
  if (include_normal) {
    out$normal <- ci_theta4_normal(
      x,
      level = level,
      se_method = se_method,
      assume_symmetric = assume_symmetric
    )
  }
  if (include_percentile) {
    out$percentile <- ci_theta4_percentile(x, B = B, level = level)
  }
  if (include_studentized) {
    out$studentized <- ci_theta4_studentized(
      x,
      B = B,
      level = level,
      se_method = se_method,
      assume_symmetric = assume_symmetric
    )
  }
  out
}

coverage_indicator <- function(ci, truth) {
  if (!is_valid_ci(ci) || !is.finite(truth)) return(NA_real_)
  as.numeric(ci$lower <= truth && truth <= ci$upper)
}

flatten_ci_result <- function(ci, truth = NA_real_) {
  if (is.null(ci)) {
    return(data.frame(
      estimate = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      width = NA_real_,
      success = FALSE,
      failure_reason = "null_ci",
      covered = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  lower <- if (!is.null(ci$lower)) ci$lower else NA_real_
  upper <- if (!is.null(ci$upper)) ci$upper else NA_real_

  data.frame(
    estimate = if (!is.null(ci$estimate)) ci$estimate else NA_real_,
    se = if (!is.null(ci$se)) ci$se else NA_real_,
    lower = lower,
    upper = upper,
    width = if (is.finite(lower) && is.finite(upper)) upper - lower else NA_real_,
    success = if (!is.null(ci$success)) isTRUE(ci$success) else FALSE,
    failure_reason = if (!is.null(ci$failure_reason)) as.character(ci$failure_reason) else NA_character_,
    covered = coverage_indicator(ci, truth),
    stringsAsFactors = FALSE
  )
}

run_one_rep_sigma2 <- function(
  n,
  dist_name,
  B = 1999L,
  level = 0.95,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE,
  se_method = c("pairwise_plugin", "empirical_projection"),
  rep_id = NA_integer_
) {
  se_method <- match.arg(se_method)

  x <- simulate_one_sample(
    n = n,
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  truth <- truth_sigma2(
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  cis <- compute_sigma2_intervals(
    x = x,
    level = level,
    B = B,
    se_method = se_method
  )

  out <- rbind(
    cbind(method = "normal", flatten_ci_result(cis$normal, truth = truth)),
    cbind(method = "percentile", flatten_ci_result(cis$percentile, truth = truth)),
    cbind(method = "studentized", flatten_ci_result(cis$studentized, truth = truth))
  )

  out$n <- n
  out$dist <- dist_name
  out$target <- "sigma2"
  out$truth <- truth
  out$rep <- rep_id
  out$se_method <- se_method
  rownames(out) <- NULL
  out
}

run_one_rep_theta4 <- function(
  n,
  dist_name,
  B = 1999L,
  level = 0.95,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE,
  se_method = c("master_formula", "empirical_projection"),
  assume_symmetric = FALSE,
  rep_id = NA_integer_
) {
  se_method <- match.arg(se_method)

  x <- simulate_one_sample(
    n = n,
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  truth <- truth_theta4(
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  cis <- compute_theta4_intervals(
    x = x,
    level = level,
    B = B,
    se_method = se_method,
    assume_symmetric = assume_symmetric
  )

  out <- rbind(
    cbind(method = "normal", flatten_ci_result(cis$normal, truth = truth)),
    cbind(method = "percentile", flatten_ci_result(cis$percentile, truth = truth)),
    cbind(method = "studentized", flatten_ci_result(cis$studentized, truth = truth))
  )

  out$n <- n
  out$dist <- dist_name
  out$target <- "theta4"
  out$truth <- truth
  out$rep <- rep_id
  out$se_method <- se_method
  out$assume_symmetric <- assume_symmetric
  rownames(out) <- NULL
  out
}

# ==========================
# Validation / debug helpers
# ==========================

check_spdv_identity <- function(x, tol = 1e-10) {
  .check_numeric_vector(x, min_n = 2L, finite_only = TRUE)
  spdv_val <- spdv(x)
  s2_val <- stats::var(x)
  list(
    spdv = spdv_val,
    var_nminus1 = s2_val,
    abs_diff = abs(spdv_val - s2_val),
    equal_within_tol = isTRUE(all.equal(spdv_val, s2_val, tolerance = tol))
  )
}

check_pairwise_estimators <- function(x) {
  .check_numeric_vector(x, min_n = 3L, finite_only = TRUE)
  list(
    n = length(x),
    sample_mean = mean(x),
    spdv = spdv(x),
    spfd = spfd(x),
    mu4_pairwise = mu4_pairwise(x),
    avar_spdv_pairwise = avar_spdv_pairwise(x),
    mu4_classical = mu4_classical(x),
    kurtosis_raw_classical = kurtosis_raw_classical(x),
    kurtosis_excess_classical = kurtosis_excess_classical(x),
    projection_spdv_mean = mean(projection_spdv(x)),
    projection_spfd_plugin_mean = mean(projection_spfd_plugin(x)),
    se_spdv_asymptotic = se_spdv_asymptotic(x),
    se_spdv_projection = se_spdv_projection(x),
    tau4_sq_plugin = tau4_sq_plugin(x),
    se_spfd_projection = se_spfd_projection(x)
  )
}

summarize_boot_failures <- function(ci_list) {
  if (is.null(ci_list) || length(ci_list) == 0L) return(data.frame())
  methods <- names(ci_list)
  data.frame(
    method = methods,
    success = vapply(ci_list, function(z) isTRUE(z$success), logical(1L)),
    failure_reason = vapply(ci_list, function(z) {
      if (is.null(z$failure_reason)) NA_character_ else as.character(z$failure_reason)
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}

# ====================
# Example master wrapper
# ====================

run_one_design <- function(
  n,
  dist_name,
  B = 1999L,
  level = 0.95,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE,
  sigma2_se_method = c("pairwise_plugin", "empirical_projection"),
  theta4_se_method = c("master_formula", "empirical_projection"),
  assume_symmetric = FALSE,
  rep_id = NA_integer_
) {
  sigma2_se_method <- match.arg(sigma2_se_method)
  theta4_se_method <- match.arg(theta4_se_method)

  x <- simulate_one_sample(
    n = n,
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  sigma_truth <- truth_sigma2(
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  theta_truth <- truth_theta4(
    dist_name = dist_name,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  sigma_ci <- compute_sigma2_intervals(
    x = x, level = level, B = B, se_method = sigma2_se_method
  )
  theta_ci <- compute_theta4_intervals(
    x = x,
    level = level,
    B = B,
    se_method = theta4_se_method,
    assume_symmetric = assume_symmetric
  )

  list(
    rep = rep_id,
    x = x,
    estimators = check_pairwise_estimators(x),
    sigma2_truth = sigma_truth,
    theta4_truth = theta_truth,
    sigma2_ci = sigma_ci,
    theta4_ci = theta_ci,
    sigma2_coverage = lapply(sigma_ci, coverage_indicator, truth = sigma_truth),
    theta4_coverage = lapply(theta_ci, coverage_indicator, truth = theta_truth)
  )
}