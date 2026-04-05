# run_all.R
# ------------------------------------------------------------------------------
# Master script to reproduce the main simulation outputs and figures
# for the pairwise U-statistics manuscript repository.
# ------------------------------------------------------------------------------

rm(list = ls())

cat("Starting full reproduction pipeline...\n")

# Create results directory if needed
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Source core helper file first
source("functions_pairwise.R")

# -----------------
# Main simulation grids
# -----------------
cat("Running sim_sigma2_grid.R ...\n")
source("sim_sigma2_grid.R")

cat("Running sim_theta4_grid.R ...\n")
source("sim_theta4_grid.R")

# -----------------
# Main figures
# -----------------
cat("Running fig_coverage_vs_df.R ...\n")
source("fig_coverage_vs_df.R")

cat("Running fig_coverage_vs_df_theta4.R ...\n")
source("fig_coverage_vs_df_theta4.R")

cat("Running fig_projection_dominance.R ...\n")
source("fig_projection_dominance.R")

# -----------------
# Session info
# -----------------
cat("Saving session information...\n")
sink("results/sessionInfo.txt")
sessionInfo()
sink()

cat("All scripts completed.\n")