# ==============================================================================
# sim_sigma2_grid.R
# Full Monte Carlo grid for sigma^2 coverage using pairwise SPDV-based intervals.
#
# Design aligned with manuscript / supplement:
#   - Target: sigma^2 via SPDV_n
#   - Methods: normal-theory, percentile bootstrap, studentized bootstrap
#   - Distributions: standardized t-laws, standard normal, contaminated normal
#   - Sample sizes: n = 50, 100, 200, 500
#   - Monte Carlo reps: M = 5000
#   - Bootstrap reps:   B = 1999
#
# Optimized for RStudio with cross-platform parallelization.
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
  library(pbapply)
})

# -----------------
# Load core methods
# -----------------
source("functions_pairwise.R")

# -----------------
# Reproducibility
# -----------------
MASTER_SEED <- 20260329L
set.seed(MASTER_SEED)

# -----------------
# Main design
# -----------------
n_grid <- c(50L, 100L, 200L, 500L)

dist_grid <- c(
  "normal", "t20", "t12", "t10", "t8", "t6", "t5", "t4", "t3", "contam"
)

M <- 5000L
B <- 1999L
level <- 0.95

se_method <- "pairwise_plugin"

eps <- 0.10
sd_clean <- 1
sd_contam <- 10
standardize_contam <- TRUE

# -----------------
# Output directories
# -----------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "sigma2_grid"), showWarnings = FALSE, recursive = TRUE)

# -----------------
# Metadata helpers
# -----------------
design_metadata <- data.frame(
  setting = c("target", "methods", "n_grid", "dist_grid", "M", "B", "level", 
              "se_method", "eps", "sd_clean", "sd_contam", "standardize_contam", "master_seed"),
  value = c("sigma2", "normal, percentile, studentized", paste(n_grid, collapse = ","), 
            paste(dist_grid, collapse = ","), as.character(M), as.character(B), 
            as.character(level), se_method, as.character(eps), as.character(sd_clean), 
            as.character(sd_contam), as.character(standardize_contam), as.character(MASTER_SEED)),
  stringsAsFactors = FALSE
)

write.csv(design_metadata, file = file.path("results", "sigma2_grid", "design_metadata_sigma2.csv"), row.names = FALSE)

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
# Design utilities
# -----------------

run_design_sigma2 <- function(
  n, dist_name, M, B, level = 0.95, eps = 0.10, sd_clean = 1, sd_contam = 10,
  standardize_contam = TRUE, se_method = "pairwise_plugin", save_each_design = TRUE,
  out_dir = file.path("results", "sigma2_grid"), cluster = NULL
) {
  
  # Export the current configuration to the worker nodes
  clusterExport(cluster, varlist = c("n", "dist_name", "B", "level", "eps", 
                                     "sd_clean", "sd_contam", "standardize_contam", 
                                     "se_method"), envir = environment())
  
  # Replace sequential loop with a parallel apply featuring an RStudio progress bar
  out_list <- pblapply(seq_len(M), function(m) {
    run_one_rep_sigma2(
      n = n, dist_name = dist_name, B = B, level = level, eps = eps,
      sd_clean = sd_clean, sd_contam = sd_contam, standardize_contam = standardize_contam,
      se_method = se_method, rep_id = m
    )
  }, cl = cluster)

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  if (save_each_design) {
    file_tag <- sprintf("sigma2_raw_%s_n%d.csv", dist_name, n)
    write.csv(out, file = file.path(out_dir, file_tag), row.names = FALSE)
  }

  out
}

summarize_sigma2_design <- function(df_design) {
  truth_val <- unique(df_design$truth)
  if (length(truth_val) != 1L) truth_val <- truth_val[1L]

  agg <- aggregate(
    cbind(covered, success, width, estimate, se, lower, upper) ~ dist + n + method + target + truth + se_method,
    data = df_design, FUN = function(z) mean(z, na.rm = TRUE)
  )

  names(agg)[names(agg) == "covered"] <- "coverage"
  names(agg)[names(agg) == "success"] <- "success_rate"
  names(agg)[names(agg) == "width"] <- "mean_width"
  names(agg)[names(agg) == "estimate"] <- "mean_estimate"
  names(agg)[names(agg) == "se"] <- "mean_se"
  names(agg)[names(agg) == "lower"] <- "mean_lower"
  names(agg)[names(agg) == "upper"] <- "mean_upper"

  fail <- aggregate(I(!success) ~ dist + n + method + target + truth + se_method, data = df_design, FUN = mean)
  names(fail)[names(fail) == "I(!success)"] <- "failure_rate"

  agg <- merge(agg, fail, by = c("dist", "n", "method", "target", "truth", "se_method"), all = TRUE)

  n_valid <- aggregate(is.finite(covered) ~ dist + n + method + target + truth + se_method, data = df_design, FUN = sum)
  names(n_valid)[names(n_valid) == "is.finite(covered)"] <- "n_coverage_nonmissing"

  agg <- merge(agg, n_valid, by = c("dist", "n", "method", "target", "truth", "se_method"), all = TRUE)

  agg$mc_reps <- aggregate(rep ~ dist + n + method + target + truth + se_method, data = df_design, FUN = length)$rep
  agg$bias <- agg$mean_estimate - agg$truth
  
  agg <- agg[order(agg$n, agg$dist, agg$method), ]
  rownames(agg) <- NULL
  agg
}

make_sigma2_table_n100 <- function(summary_df) {
  tab <- subset(summary_df, n == 100L, select = c("dist", "method", "coverage"))
  wide <- reshape(tab, idvar = "dist", timevar = "method", direction = "wide")
  
  names(wide) <- sub("^coverage\\.", "", names(wide))
  names(wide) <- sub("^method\\.", "", names(wide))
  wide <- wide[order(match(wide$dist, dist_grid)), ]
  rownames(wide) <- NULL
  wide
}

# -----------------
# Run full grid
# -----------------
all_results_list <- vector("list", length(n_grid) * length(dist_grid))
all_summary_list <- vector("list", length(n_grid) * length(dist_grid))

k <- 1L

for (n in n_grid) {
  for (dist_name in dist_grid) {

    message("--------------------------------------------------")
    message(sprintf("Starting sigma2 design: dist=%s, n=%d", dist_name, n))
    
    design_seed <- MASTER_SEED + 100000L * match(as.character(n), as.character(n_grid)) +
      1000L * match(dist_name, dist_grid)
    set.seed(design_seed)

    # Pass the pre-initialized cluster into the grid runner
    df_design <- run_design_sigma2(
      n = n, dist_name = dist_name, M = M, B = B, level = level, eps = eps,
      sd_clean = sd_clean, sd_contam = sd_contam, standardize_contam = standardize_contam,
      se_method = se_method, save_each_design = TRUE, out_dir = file.path("results", "sigma2_grid"),
      cluster = cl
    )

    summary_design <- summarize_sigma2_design(df_design)
    write.csv(summary_design, file = file.path("results", "sigma2_grid", sprintf("sigma2_summary_%s_n%d.csv", dist_name, n)), row.names = FALSE)

    all_results_list[[k]] <- df_design
    all_summary_list[[k]] <- summary_design
    k <- k + 1L
  }
}

# Safely tear down cluster when finished
stopCluster(cl)

res_sigma2_all <- do.call(rbind, all_results_list)
summary_sigma2_all <- do.call(rbind, all_summary_list)

rownames(res_sigma2_all) <- NULL
rownames(summary_sigma2_all) <- NULL

# -----------------
# Save combined files
# -----------------
write.csv(res_sigma2_all, file = file.path("results", "sigma2_grid", "sigma2_all_raw.csv"), row.names = FALSE)
write.csv(summary_sigma2_all, file = file.path("results", "sigma2_grid", "sigma2_all_summary.csv"), row.names = FALSE)

table_sigma2_n100 <- make_sigma2_table_n100(summary_sigma2_all)
write.csv(table_sigma2_n100, file = file.path("results", "sigma2_grid", "sigma2_table_n100.csv"), row.names = FALSE)

summary_sigma2_compact <- subset(
  summary_sigma2_all,
  select = c("dist", "n", "method", "truth", "coverage", "failure_rate", "success_rate",
             "mean_width", "mean_estimate", "bias", "mean_se", "mc_reps", "n_coverage_nonmissing", "se_method")
)
write.csv(summary_sigma2_compact, file = file.path("results", "sigma2_grid", "sigma2_all_summary_compact.csv"), row.names = FALSE)

# -----------------
# Console output
# -----------------
message("==================================================")
message("Full sigma2 simulation complete.")
print(table_sigma2_n100)