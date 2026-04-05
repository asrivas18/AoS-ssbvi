# fig_coverage_vs_df.R
# ------------------------------------------------------------------------------
# Coverage curves versus degrees of freedom at fixed n = 100
#
# Supplement-style figure:
# Coverage curves as a function of df for standardized t distributions,
# showing threshold behavior near df = 4, 8, 12 for 95% confidence intervals.
#
# This script focuses on sigma^2 inference using SPDV_n, comparing:
#   - Normal-theory interval
#   - Percentile bootstrap interval
#   - Studentized bootstrap interval
#
# Uses the snake_case helper API from functions_pairwise.R.
#
# Outputs:
# - results/fig_coverage_vs_df_raw.csv
# - results/fig_coverage_vs_df_summary.csv
# - results/fig_coverage_vs_df.pdf
# - results/fig_coverage_vs_df.png
# ------------------------------------------------------------------------------

rm(list = ls())

source("functions_pairwise.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20260330)

# -----------------
# Design
# -----------------
n <- 100L
M <- 1000L          # raise to 5000L for publication-quality final runs
B <- 499L           # raise to 1999L for publication-quality final runs
level <- 0.95

# Dense df grid for threshold plot
df_grid <- c(3, 4, 5, 6, 8, 10, 12, 20)

dist_grid <- paste0("t", df_grid)

# -----------------
# One-rep wrapper
# -----------------
one_rep_sigma2 <- function(n, dist_name, B, level, rep_id) {
  out <- run_one_rep_sigma2(
    n = n,
    dist_name = dist_name,
    B = B,
    level = level,
    se_method = "pairwise_plugin",
    rep_id = rep_id
  )

  out$df <- as.numeric(sub("^t", "", dist_name))
  out
}

# -----------------
# Simulation
# -----------------
out_list <- vector("list", length(dist_grid) * M)
k <- 1L

for (dist_name in dist_grid) {
  message(sprintf("Running sigma^2 coverage: dist=%s, n=%d, M=%d, B=%d",
                  dist_name, n, M, B))

  for (m in seq_len(M)) {
    if (m %% 100L == 0L) {
      message(sprintf("  rep %d / %d", m, M))
    }

    out_list[[k]] <- one_rep_sigma2(
      n = n,
      dist_name = dist_name,
      B = B,
      level = level,
      rep_id = m
    )
    k <- k + 1L
  }
}

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
# Coverage summary
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
    p <- mean(z, na.rm = TRUE)
    m <- sum(is.finite(z))
    sqrt(p * (1 - p) / m)
  }
)

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
  geom_hline(
    yintercept = level,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey35"
  ) +
  geom_vline(
    xintercept = c(4, 8, 12),
    linetype = "dotted",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_ribbon(
    aes(ymin = lower_mc, ymax = upper_mc, fill = method),
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = df_grid) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
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
  annotate("text", x = 4, y = 0.05, label = "df = 4", angle = 90, vjust = -0.4, size = 3.3) +
  annotate("text", x = 8, y = 0.05, label = "df = 8", angle = 90, vjust = -0.4, size = 3.3) +
  annotate("text", x = 12, y = 0.05, label = "df = 12", angle = 90, vjust = -0.4, size = 3.3) +
  labs(
    title = expression(paste("Coverage vs Degrees of Freedom for ", sigma^2, " Intervals")),
    subtitle = paste0(
      "Standardized t distributions, n = ", n,
      ", Monte Carlo replications = ", M,
      ", bootstrap resamples = ", B
    ),
    x = "Degrees of freedom",
    y = "Empirical coverage",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "results/fig_coverage_vs_df.pdf",
  plot = p,
  width = 9.2,
  height = 6.6
)

ggsave(
  filename = "results/fig_coverage_vs_df.png",
  plot = p,
  width = 9.2,
  height = 6.6,
  dpi = 320
)

# -----------------
# Save outputs
# -----------------
write.csv(
  res,
  file = "results/fig_coverage_vs_df_raw.csv",
  row.names = FALSE
)

write.csv(
  summary_tab,
  file = "results/fig_coverage_vs_df_summary.csv",
  row.names = FALSE
)

print(summary_tab)
print(p)