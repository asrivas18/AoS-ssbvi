# ==============================================================================
# fig_coverage_transition_sigma2.R
# Coverage curves versus degrees of freedom for sigma^2 (Figure 4)
#
# Empirical coverage of nominal 95% confidence intervals for sigma^2 based on SPDV_n
# across standardized t distributions, highlighting the threshold at df = 4.
# Matches supp_coverage_transition_sigma2.pdf in the LaTeX supplement.
# ==============================================================================

rm(list = ls())

source("functions_pairwise.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(parallel)
  library(pbapply) # For parallel progress bars
})

set.seed(20260505)

# -----------------
# Design
# -----------------
n <- 100L
M <- 5000L          # High replication for stable MCSE bands
B <- 1999L          # High bootstrap replication for studentized stability
level <- 0.95

# Supplement-style threshold grid
df_grid <- c(3, 4, 5, 6, 8, 10, 12, 20)
dist_grid <- paste0("t", df_grid)

# Default studentization choice for sigma^2
se_method <- "pairwise_plugin"

# -----------------
# Setup Cross-Platform Parallel Cluster
# -----------------
cores <- max(1L, detectCores() - 1L)
message(sprintf("Setting up parallel cluster with %d cores...", cores))

cl <- makeCluster(cores)
clusterEvalQ(cl, {
  source("functions_pairwise.R")
})

# -----------------
# Simulation Loop
# -----------------
out_list <- vector("list", length(dist_grid))
names(out_list) <- dist_grid

for (i in seq_along(dist_grid)) {
  dist_name <- dist_grid[i]
  message(sprintf("Running sigma2 coverage: dist=%s, n=%d, M=%d, B=%d", dist_name, n, M, B))
  
  # Export variables to the cluster for this specific distribution
  clusterExport(cl, varlist = c("n", "dist_name", "B", "level", "se_method"), envir = environment())
  
  # Run M replications in parallel
  reps_list <- pblapply(seq_len(M), function(m) {
    out <- run_one_rep_sigma2(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      se_method = se_method,
      rep_id = m
    )
    out$df <- as.numeric(sub("^t", "", dist_name))
    out
  }, cl = cl)
  
  out_list[[i]] <- do.call(rbind, reps_list)
}

# Safely shut down cluster
stopCluster(cl)

res <- do.call(rbind, out_list)
rownames(res) <- NULL

# -----------------
# Clean labels
# -----------------
res$method <- factor(
  res$method,
  levels = c("normal", "percentile", "studentized"),
  labels = c("Normal-theory", "Percentile bootstrap", "Studentized bootstrap")
)

# -----------------
# Coverage summary & MCSE
# -----------------
coverage_tab <- aggregate(covered ~ df + method, data = res, FUN = function(z) mean(z, na.rm = TRUE))
names(coverage_tab)[names(coverage_tab) == "covered"] <- "coverage"

mcse_tab <- aggregate(covered ~ df + method, data = res, FUN = function(z) {
  zf <- z[is.finite(z)]
  p <- mean(zf)
  m <- length(zf)
  if(m == 0) return(NA_real_)
  sqrt(p * (1 - p) / m)
})
names(mcse_tab)[names(mcse_tab) == "covered"] <- "mcse"

summary_tab <- merge(coverage_tab, mcse_tab, by = c("df", "method"), sort = TRUE)

summary_tab$nominal <- level
summary_tab$lower_mc <- pmax(summary_tab$coverage - 1.96 * summary_tab$mcse, 0)
summary_tab$upper_mc <- pmin(summary_tab$coverage + 1.96 * summary_tab$mcse, 1)

# -----------------
# Plot
# -----------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

p <- ggplot(summary_tab, aes(x = df, y = coverage, color = method)) +
  geom_hline(yintercept = level, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  
  # The critical 4th moment boundary for SPDV_n
  geom_vline(xintercept = 4, linetype = "dotted", linewidth = 0.6, color = "grey40") +
  
  geom_ribbon(aes(ymin = lower_mc, ymax = upper_mc, fill = method), alpha = 0.14, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.5) +
  
  # Log10 scaling stretches the critical transition zone
  scale_x_log10(breaks = df_grid, minor_breaks = NULL) +
  
  # y-axis set to capture the specific breakdown severity of sigma^2 inference
  scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, by = 0.1)) +
  
  scale_color_manual(
    values = c("Normal-theory" = "#4C78A8", "Percentile bootstrap" = "#54A24B", "Studentized bootstrap" = "#E45756")
  ) +
  scale_fill_manual(
    values = c("Normal-theory" = "#4C78A8", "Percentile bootstrap" = "#54A24B", "Studentized bootstrap" = "#E45756")
  ) +
  
  # Annotation for the boundary
  annotate("text", x = 4.2, y = 0.45, label = "4th moment threshold", angle = 90, vjust = 1, size = 4.2, color="grey30") +
  
  labs(
    title = expression(paste("Coverage vs Degrees of Freedom for ", sigma^2, " Intervals")),
    subtitle = paste0("Inference based on ", SPDV[n], " under standardized t distributions (n = ", n, ")"),
    x = "Degrees of freedom (log scale)",
    y = "Empirical coverage probability",
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
    axis.text = element_text(color = "black")
  )

# -----------------
# Save outputs matching the LaTeX \includegraphics exactly
# -----------------
ggsave(filename = "results/supp_coverage_transition_sigma2.pdf", plot = p, width = 9.2, height = 6.6)
ggsave(filename = "results/supp_coverage_transition_sigma2.png", plot = p, width = 9.2, height = 6.6, dpi = 320)

write.csv(res, file = "results/supp_coverage_transition_sigma2_raw.csv", row.names = FALSE)
write.csv(summary_tab, file = "results/supp_coverage_transition_sigma2_summary.csv", row.names = FALSE)

cat("\n=== Sigma^2 Coverage Transition Complete ===\n")
print(summary_tab[, c("df", "method", "coverage", "mcse")])
cat("Files saved in 'results/' folder.\n")