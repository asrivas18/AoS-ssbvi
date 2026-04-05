# sim_theta4_full.R
# ------------------------------------------------------------------------------
# Full Monte Carlo grid for theta4 coverage using SPFD-based intervals.
#
# Matches the manuscript/supplement design:
# - n in {50, 100, 200, 500}
# - distributions: normal, t3, t4, t5, t6, t8, t10, t12, t20, contam
# - M = 5000 Monte Carlo replications
# - B = 1999 bootstrap resamples per dataset
#
# Outputs:
# - results/res_theta4_full.csv
# - results/summary_theta4_full.csv
# ------------------------------------------------------------------------------

rm(list = ls())

source("functions_pairwise.R")

set.seed(20260330)

# -----------------
# Full-grid design
# -----------------
n_grid <- c(50L, 100L, 200L, 500L)
dist_grid <- c("normal", "t3", "t4", "t5", "t6", "t8", "t10", "t12", "t20", "contam")

M <- 5000L
B <- 1999L
level <- 0.95

# contaminated-normal settings
eps <- 0.10
sd_clean <- 1
sd_contam <- 10
standardize_contam <- TRUE

# theta4 interval settings
se_method <- "master_formula"
assume_symmetric <- TRUE

# -----------------
# One-design runner
# -----------------
run_design_theta4_full <- function(
  n,
  dist_name,
  M,
  B,
  level = 0.95,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE,
  se_method = "master_formula",
  assume_symmetric = TRUE
) {
  out_list <- vector("list", M)

  for (m in seq_len(M)) {
    if (m %% 100 == 0L) {
      message(sprintf("  rep %d / %d", m, M))
    }

    out_list[[m]] <- run_one_rep_theta4(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam,
      se_method = se_method,
      assume_symmetric = assume_symmetric,
      rep_id = m
    )
  }

  do.call(rbind, out_list)
}

# -------------------------
# Run full simulation grid
# -------------------------
res_list <- vector("list", length(n_grid) * length(dist_grid))
k <- 1L

for (n in n_grid) {
  for (dist_name in dist_grid) {
    message(sprintf(
      "Running theta4 full grid: n=%d, dist=%s, M=%d, B=%d",
      n, dist_name, M, B
    ))

    res_list[[k]] <- run_design_theta4_full(
      n = n,
      dist_name = dist_name,
      M = M,
      B = B,
      level = level,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam,
      se_method = se_method,
      assume_symmetric = assume_symmetric
    )

    k <- k + 1L
  }
}

res_theta4_full <- do.call(rbind, res_list)
rownames(res_theta4_full) <- NULL

# -------------------------
# Summaries
# -------------------------
res_theta4_full$covered_uncond <- ifelse(
  is.na(res_theta4_full$covered), 0, res_theta4_full$covered
)

res_theta4_full$covered_cond <- with(
  res_theta4_full,
  ifelse(success, covered, NA_real_)
)

res_theta4_full$valid_ci <- with(
  res_theta4_full,
  is.finite(lower) & is.finite(upper) & (lower <= upper)
)

summary_theta4_full <- aggregate(
  cbind(
    covered_uncond,
    covered_cond,
    success,
    valid_ci,
    width
  ) ~ dist + n + method + target + se_method + assume_symmetric,
  data = res_theta4_full,
  FUN = function(z) mean(z, na.rm = TRUE)
)

names(summary_theta4_full)[names(summary_theta4_full) == "covered_uncond"] <- "coverage_uncond"
names(summary_theta4_full)[names(summary_theta4_full) == "covered_cond"] <- "coverage_cond"
names(summary_theta4_full)[names(summary_theta4_full) == "success"] <- "success_rate"
names(summary_theta4_full)[names(summary_theta4_full) == "valid_ci"] <- "valid_ci_rate"
names(summary_theta4_full)[names(summary_theta4_full) == "width"] <- "mean_width"

failure_theta4_full <- aggregate(
  I(!success) ~ dist + n + method + target + se_method + assume_symmetric,
  data = res_theta4_full,
  FUN = mean
)

names(failure_theta4_full)[names(failure_theta4_full) == "I(!success)"] <- "failure_rate"

summary_theta4_full <- merge(
  summary_theta4_full,
  failure_theta4_full,
  by = c("dist", "n", "method", "target", "se_method", "assume_symmetric"),
  all = TRUE
)

summary_theta4_full <- summary_theta4_full[
  order(summary_theta4_full$n, summary_theta4_full$dist, summary_theta4_full$method),
]

rownames(summary_theta4_full) <- NULL

# -------------------------
# Save outputs
# -------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(
  res_theta4_full,
  file = "results/res_theta4_full.csv",
  row.names = FALSE
)

write.csv(
  summary_theta4_full,
  file = "results/summary_theta4_full.csv",
  row.names = FALSE
)

print(summary_theta4_full)