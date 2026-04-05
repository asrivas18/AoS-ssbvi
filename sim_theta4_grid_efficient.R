# ==============================================================================
# sim_theta4_grid.R
# Full Monte Carlo grid for theta4 coverage using SPFD-based intervals.
#
# Optimized for RStudio with cross-platform parallelization.
#
# Matches the manuscript/supplement design:
# - n in {50, 100, 200, 500}
# - distributions: normal, t3, t4, t5, t6, t8, t10, t12, t20, contam
# - M = 5000 Monte Carlo replications
# - B = 1999 bootstrap resamples per dataset
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
  library(pbapply)
})

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

# -----------------
# Setup Parallel Cluster
# -----------------
cores <- max(1L, detectCores() - 1L)
message(sprintf("Setting up parallel cluster with %d cores...", cores))

cl <- makeCluster(cores)
clusterEvalQ(cl, {
  source("functions_pairwise.R")
})

# -----------------
# One-design runner
# -----------------
run_design_theta4_full <- function(
  n, dist_name, M, B, level = 0.95, eps = 0.10, sd_clean = 1, 
  sd_contam = 10, standardize_contam = TRUE, cluster = NULL
) {
  
  # Export configuration variables to the worker nodes
  clusterExport(cluster, varlist = c("n", "dist_name", "B", "level", "eps", 
                                     "sd_clean", "sd_contam", "standardize_contam"), 
                envir = environment())

  # Run M replications in parallel with a progress bar
  out_list <- pblapply(seq_len(M), function(m) {
    # FIX: Removed unsupported se_method and assume_symmetric arguments
    run_one_rep_theta4(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam,
      rep_id = m
    )
  }, cl = cluster)

  do.call(rbind, out_list)
}

# -------------------------
# Run full simulation grid
# -------------------------
res_list <- vector("list", length(n_grid) * length(dist_grid))
k <- 1L

for (n in n_grid) {
  for (dist_name in dist_grid) {
    message("--------------------------------------------------")
    message(sprintf("Running theta4 full grid: dist=%s, n=%d, M=%d, B=%d", dist_name, n, M, B))

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
      cluster = cl
    )

    k <- k + 1L
  }
}

# Safely shut down cluster when finished
stopCluster(cl)

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

# FIX: Removed se_method and assume_symmetric from aggregation logic
summary_theta4_full <- aggregate(
  cbind(covered_uncond, covered_cond, success, valid_ci, width) ~ dist + n + method + target,
  data = res_theta4_full,
  FUN = function(z) mean(z, na.rm = TRUE)
)

names(summary_theta4_full)[names(summary_theta4_full) == "covered_uncond"] <- "coverage_uncond"
names(summary_theta4_full)[names(summary_theta4_full) == "covered_cond"] <- "coverage_cond"
names(summary_theta4_full)[names(summary_theta4_full) == "success"] <- "success_rate"
names(summary_theta4_full)[names(summary_theta4_full) == "valid_ci"] <- "valid_ci_rate"
names(summary_theta4_full)[names(summary_theta4_full) == "width"] <- "mean_width"

failure_theta4_full <- aggregate(
  I(!success) ~ dist + n + method + target,
  data = res_theta4_full,
  FUN = mean
)

names(failure_theta4_full)[names(failure_theta4_full) == "I(!success)"] <- "failure_rate"

summary_theta4_full <- merge(
  summary_theta4_full,
  failure_theta4_full,
  by = c("dist", "n", "method", "target"),
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

message("==================================================")
message("Full theta4 simulation complete.")
print(head(summary_theta4_full))