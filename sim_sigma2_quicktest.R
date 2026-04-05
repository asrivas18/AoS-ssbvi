# sim_sigma2_quicktest.R
# ------------------------------------------------------------------------------
# Quick pilot simulation for sigma^2 coverage using functions_pairwise.R
# Intended for debugging before the full Monte Carlo grid.
# ------------------------------------------------------------------------------

rm(list = ls())

source("functions_pairwise.R")

set.seed(20260329)

# -----------------
# Quick-test design
# -----------------
n_grid <- c(50, 100)
dist_grid <- c("normal", "t4", "t8", "t12", "contam")

M <- 50L        # quick pilot; increase to 100 if desired
B <- 199L       # quick pilot bootstrap count
level <- 0.95
se_method <- "pairwise_plugin"

# contaminated normal settings used in manuscript/supplement style
eps <- 0.10
sd_clean <- 1
sd_contam <- 10
standardize_contam <- TRUE

# -----------------
# One-design runner
# -----------------
run_design_sigma2_quick <- function(
  n,
  dist_name,
  M,
  B,
  level = 0.95,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE,
  se_method = "pairwise_plugin"
) {
  out_list <- vector("list", M)

  for (m in seq_len(M)) {
    out_list[[m]] <- run_one_rep_sigma2(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam,
      se_method = se_method,
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
    message(sprintf("Running quick test: n=%d, dist=%s", n, dist_name))
    res_list[[k]] <- run_design_sigma2_quick(
      n = n,
      dist_name = dist_name,
      M = M,
      B = B,
      level = level,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam,
      se_method = se_method
    )
    k <- k + 1L
  }
}

res_sigma2_quick <- do.call(rbind, res_list)
rownames(res_sigma2_quick) <- NULL

# -------------------------
# Summaries
# -------------------------
summary_sigma2_quick <- aggregate(
  cbind(
    covered,
    success,
    width
  ) ~ dist + n + method + target + se_method,
  data = res_sigma2_quick,
  FUN = function(z) mean(z, na.rm = TRUE)
)

names(summary_sigma2_quick)[names(summary_sigma2_quick) == "covered"] <- "coverage"
names(summary_sigma2_quick)[names(summary_sigma2_quick) == "success"] <- "success_rate"
names(summary_sigma2_quick)[names(summary_sigma2_quick) == "width"] <- "mean_width"

failure_sigma2_quick <- aggregate(
  I(!success) ~ dist + n + method + target + se_method,
  data = res_sigma2_quick,
  FUN = mean
)

names(failure_sigma2_quick)[names(failure_sigma2_quick) == "I(!success)"] <- "failure_rate"

# merge
summary_sigma2_quick <- merge(
  summary_sigma2_quick,
  failure_sigma2_quick,
  by = c("dist", "n", "method", "target", "se_method"),
  all = TRUE
)

# sort for readability
summary_sigma2_quick <- summary_sigma2_quick[
  order(summary_sigma2_quick$n, summary_sigma2_quick$dist, summary_sigma2_quick$method),
]

# -------------------------
# Save outputs
# -------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(
  res_sigma2_quick,
  file = "results/res_sigma2_quick.csv",
  row.names = FALSE
)

write.csv(
  summary_sigma2_quick,
  file = "results/summary_sigma2_quick.csv",
  row.names = FALSE
)

print(summary_sigma2_quick)