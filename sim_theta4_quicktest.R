# sim_theta4_quicktest.R
# ------------------------------------------------------------------------------
# Quick pilot simulation for theta4 coverage using SPFD-based intervals.
# Intended for debugging before the full Monte Carlo grid.
# ------------------------------------------------------------------------------

rm(list = ls())

source("functions_pairwise.R")

set.seed(20260329)

# -----------------
# Quick-test design
# -----------------
n_grid <- c(50L, 100L)
dist_grid <- c("normal", "t6", "t8", "t12", "contam")

M <- 50L
B <- 199L
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
run_design_theta4_quick <- function(
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
# Run quick simulation grid
# -------------------------
res_list <- vector("list", length(n_grid) * length(dist_grid))
k <- 1L

for (n in n_grid) {
  for (dist_name in dist_grid) {
    message(sprintf("Running theta4 quick test: n=%d, dist=%s", n, dist_name))
    res_list[[k]] <- run_design_theta4_quick(
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

res_theta4_quick <- do.call(rbind, res_list)
rownames(res_theta4_quick) <- NULL

# -------------------------
# Summaries
# -------------------------
res_theta4_quick$covered_uncond <- ifelse(
  is.na(res_theta4_quick$covered), 0, res_theta4_quick$covered
)

res_theta4_quick$covered_cond <- with(
  res_theta4_quick,
  ifelse(success, covered, NA_real_)
)

res_theta4_quick$valid_ci <- with(
  res_theta4_quick,
  is.finite(lower) & is.finite(upper) & (lower <= upper)
)

summary_theta4_quick <- aggregate(
  cbind(
    covered_uncond,
    covered_cond,
    success,
    valid_ci,
    width
  ) ~ dist + n + method + target + se_method + assume_symmetric,
  data = res_theta4_quick,
  FUN = function(z) mean(z, na.rm = TRUE)
)

names(summary_theta4_quick)[names(summary_theta4_quick) == "covered_uncond"] <- "coverage_uncond"
names(summary_theta4_quick)[names(summary_theta4_quick) == "covered_cond"] <- "coverage_cond"
names(summary_theta4_quick)[names(summary_theta4_quick) == "success"] <- "success_rate"
names(summary_theta4_quick)[names(summary_theta4_quick) == "valid_ci"] <- "valid_ci_rate"
names(summary_theta4_quick)[names(summary_theta4_quick) == "width"] <- "mean_width"

failure_theta4_quick <- aggregate(
  I(!success) ~ dist + n + method + target + se_method + assume_symmetric,
  data = res_theta4_quick,
  FUN = mean
)

names(failure_theta4_quick)[names(failure_theta4_quick) == "I(!success)"] <- "failure_rate"

summary_theta4_quick <- merge(
  summary_theta4_quick,
  failure_theta4_quick,
  by = c("dist", "n", "method", "target", "se_method", "assume_symmetric"),
  all = TRUE
)

summary_theta4_quick <- summary_theta4_quick[
  order(summary_theta4_quick$n, summary_theta4_quick$dist, summary_theta4_quick$method),
]

# -------------------------
# Save outputs
# -------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(
  res_theta4_quick,
  file = "results/res_theta4_quick.csv",
  row.names = FALSE
)

write.csv(
  summary_theta4_quick,
  file = "results/summary_theta4_quick.csv",
  row.names = FALSE
)

print(summary_theta4_quick)