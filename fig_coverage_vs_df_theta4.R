# ==============================================================================
# fig_coverage_vs_df_theta4.R
# Coverage curves versus degrees of freedom at fixed n = 100
#
# Empirical coverage of nominal 95% confidence intervals for theta4 based on SPFD_n
# across standardized t distributions, highlighting thresholds at df = 4, 8, 12.
# ==============================================================================

rm(list = ls())

source("functions_pairwise.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(parallel)
  library(pbapply)
})

set.seed(20260330)

# -----------------
# Design
# -----------------
n <- 100L
M <- 1000L
B <- 499L
level <- 0.95

df_grid <- c(3, 4, 5, 6, 8, 10, 12, 20)
dist_grid <- paste0("t", df_grid)

# -----------------
# Parallel cluster
# -----------------
cores <- max(1L, detectCores() - 1L)
message(sprintf("Setting up parallel cluster with %d cores...", cores))

cl <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  source("functions_pairwise.R")
  NULL
})

clusterSetRNGStream(cl, iseed = 20260330)

# -----------------
# Simulation loop
# -----------------
out_list <- vector("list", length(dist_grid))
names(out_list) <- dist_grid

for (i in seq_along(dist_grid)) {
  dist_name <- dist_grid[i]
  message(sprintf(
    "Running theta4 coverage: dist=%s, n=%d, M=%d, B=%d",
    dist_name, n, M, B
  ))

  clusterExport(
    cl,
    varlist = c("n", "dist_name", "B", "level"),
    envir = environment()
  )

  reps_list <- pblapply(seq_len(M), function(m) {
    ans <- run_one_design(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      theta4_se_method = "master_formula",
      assume_symmetric = FALSE,
      rep_id = m
    )

    truth <- ans$theta4_truth

    out <- rbind(
      cbind(method = "normal",      flatten_ci_result(ans$theta4_ci$normal,      truth = truth)),
      cbind(method = "percentile",  flatten_ci_result(ans$theta4_ci$percentile,  truth = truth)),
      cbind(method = "studentized", flatten_ci_result(ans$theta4_ci$studentized, truth = truth))
    )

    out$n <- n
    out$dist <- dist_name
    out$target <- "theta4"
    out$truth <- truth
    out$rep <- m
    out$theta4_se_method <- "master_formula"
    out$assume_symmetric <- FALSE
    out$df <- as.numeric(sub("^t", "", dist_name))
    rownames(out) <- NULL
    out
  }, cl = cl)

  out_list[[i]] <- do.call(rbind, reps_list)
}

res <- do.call(rbind, out_list)
rownames(res) <- NULL

# -----------------
# Labels
# -----------------
res$method <- factor(
  res$method,
  levels = c("normal", "percentile", "studentized"),
  labels = c("Normal-theory", "Percentile bootstrap", "Studentized bootstrap")
)

# numeric conversions after cbind
num_cols <- c("estimate", "se", "lower", "upper", "width", "covered", "truth", "df")
for (nm in intersect(num_cols, names(res))) {
  res[[nm]] <- as.numeric(res[[nm]])
}

# -----------------
# Summary
# -----------------
coverage_tab <- aggregate(
  covered ~ df + method,
  data = res,
  FUN = function(z) mean(z, na.rm = TRUE)
)
names(coverage_tab)[names(coverage_tab) == "covered"] <- "coverage"

mcse_tab <- aggregate(
  covered ~ df + method,
  data = res,
  FUN = function(z) {
    zf <- z[is.finite(z)]
    p <- mean(zf)
    m <- length(zf)
    sqrt(p * (1 - p) / m)
  }
)
names(mcse_tab)[names(mcse_tab) == "covered"] <- "mcse"

avg_len_tab <- aggregate(
  width ~ df + method,
  data = res,
  FUN = function(z) mean(z, na.rm = TRUE)
)
names(avg_len_tab)[names(avg_len_tab) == "width"] <- "avg_width"

summary_tab <- Reduce(
  function(x, y) merge(x, y, by = c("df", "method"), sort = TRUE),
  list(coverage_tab, mcse_tab, avg_len_tab)
)

summary_tab$nominal <- level
summary_tab$lower_mc <- pmax(summary_tab$coverage - 1.96 * summary_tab$mcse, 0)
summary_tab$upper_mc <- pmin(summary_tab$coverage + 1.96 * summary_tab$mcse, 1)

# -----------------
# Plot
# -----------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

p <- ggplot(summary_tab, aes(x = df, y = coverage, color = method)) +
  geom_hline(
    yintercept = level,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey35"
  ) +
  geom_vline(
    xintercept = c(4, 8, 12),
    linetype = "dotted",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_ribbon(
    aes(ymin = lower_mc, ymax = upper_mc, fill = method),
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = df_grid) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(
    values = c(
      "Normal-theory" = "#4C78A8",
      "Percentile bootstrap" = "#54A24B",
      "Studentized bootstrap" = "#E45756"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Normal-theory" = "#4C78A8",
      "Percentile bootstrap" = "#54A24B",
      "Studentized bootstrap" = "#E45756"
    )
  ) +
  annotate("text", x = 4, y = 0.05, label = "4th moment", angle = 90, vjust = -0.4, size = 3.6, color = "grey30") +
  annotate("text", x = 8, y = 0.05, label = "8th moment", angle = 90, vjust = -0.4, size = 3.6, color = "grey30") +
  annotate("text", x = 12, y = 0.05, label = "12th moment", angle = 90, vjust = -0.4, size = 3.6, color = "grey30") +
  labs(
    title = expression(paste("Coverage vs Degrees of Freedom for ", theta[4], " Intervals")),
    subtitle = paste0(
      "Standardized t distributions, n = ", n,
      ", MC replications = ", M,
      ", bootstrap resamples = ", B
    ),
    x = "Degrees of freedom",
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

ggsave("results/fig_coverage_vs_df_theta4.pdf", plot = p, width = 9.2, height = 6.6)
ggsave("results/fig_coverage_vs_df_theta4.png", plot = p, width = 9.2, height = 6.6, dpi = 320)

write.csv(res, "results/fig_coverage_vs_df_theta4_raw.csv", row.names = FALSE)
write.csv(summary_tab, "results/fig_coverage_vs_df_theta4_summary.csv", row.names = FALSE)

print(summary_tab[, c("df", "method", "coverage", "mcse", "avg_width")])