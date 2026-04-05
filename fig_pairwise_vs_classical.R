# fig_pairwise_vs_classical.R
# ------------------------------------------------------------------------------
# Publication-ready ggplot2 figure:
# Pairwise vs classical estimation of mu4
#
# Figure S1 style:
# Boxplots comparing
#   - Pairwise estimator: 6*SPFD_n - 3*SPDV_n^2 = mu4_pairwise(x)
#   - Classical estimator: mean((x - mean(x))^4) = mu4_classical(x)
# across selected standardized t distributions and contaminated normal.
#
# Outputs:
# - results/fig_pairwise_vs_classical_raw.csv
# - results/fig_pairwise_vs_classical_summary.csv
# - results/fig_pairwise_vs_classical.pdf
# - results/fig_pairwise_vs_classical.png
# ------------------------------------------------------------------------------

rm(list = ls())

source("functions_pairwise.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20260330)

# -----------------
# Figure design
# -----------------
n <- 100L
M <- 2000L

# Selected designs referenced in supplement-style estimator comparison
dist_grid <- c("t20", "t10", "t6", "contam")

dist_labels <- c(
  t20 = "t(20)",
  t10 = "t(10)",
  t6 = "t(6)",
  contam = "Contaminated normal"
)

# contaminated-normal settings
eps <- 0.10
sd_clean <- 1
sd_contam <- 10
standardize_contam <- TRUE

# -----------------
# One-rep generator
# -----------------
one_rep_pairwise_vs_classical <- function(
  n,
  distname,
  eps = 0.10,
  sd_clean = 1,
  sd_contam = 10,
  standardize_contam = TRUE
) {
  # FIX: Corrected to match exact function names in functions_pairwise.R
  x <- simulate_one_sample(
    n = n,
    dist_name = distname,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  mu4_true <- truth_mu4(
    dist_name = distname,
    eps = eps,
    sd_clean = sd_clean,
    sd_contam = sd_contam,
    standardize_contam = standardize_contam
  )

  est_pairwise <- mu4_pairwise(x)
  est_classical <- mu4_classical(x)

  data.frame(
    dist = distname,
    n = n,
    truth_mu4 = mu4_true,
    estimator = c("Pairwise", "Classical"),
    estimate = c(est_pairwise, est_classical),
    bias = c(est_pairwise - mu4_true, est_classical - mu4_true),
    sq_error = c((est_pairwise - mu4_true)^2, (est_classical - mu4_true)^2),
    stringsAsFactors = FALSE
  )
}

# -----------------
# Run simulation
# -----------------
out_list <- vector("list", length(dist_grid) * M)
k <- 1L

for (distname in dist_grid) {
  message(sprintf("Running estimator comparison: dist=%s, n=%d, M=%d", distname, n, M))

  for (m in seq_len(M)) {
    if (m %% 200L == 0L) {
      message(sprintf("  rep %d / %d", m, M))
    }

    out_list[[k]] <- one_rep_pairwise_vs_classical(
      n = n,
      distname = distname,
      eps = eps,
      sd_clean = sd_clean,
      sd_contam = sd_contam,
      standardize_contam = standardize_contam
    )
    k <- k + 1L
  }
}

res <- do.call(rbind, out_list)
rownames(res) <- NULL

# -----------------
# Factor cleanup
# -----------------
res$dist <- factor(res$dist, levels = dist_grid, labels = unname(dist_labels[dist_grid]))
res$estimator <- factor(res$estimator, levels = c("Pairwise", "Classical"))

# -----------------
# Summary table
# -----------------
summ_fun <- function(z) {
  c(
    mean = mean(z, na.rm = TRUE),
    sd = stats::sd(z, na.rm = TRUE),
    median = stats::median(z, na.rm = TRUE)
  )
}

summary_tab <- aggregate(
  cbind(estimate, bias, sq_error) ~ dist + n + estimator + truth_mu4,
  data = res,
  FUN = summ_fun
)

summary_tab <- do.call(data.frame, summary_tab)
names(summary_tab) <- gsub("\\.", "_", names(summary_tab))

if ("sq_error_mean" %in% names(summary_tab)) {
  summary_tab$rmse <- sqrt(summary_tab$sq_error_mean)
}

# facet-level truth labels
truth_df <- unique(res[, c("dist", "truth_mu4")])
truth_df$truth_lab <- paste0("Truth = ", format(round(truth_df$truth_mu4, 3), nsmall = 3))

# -----------------
# Plot
# -----------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 
p <- ggplot(res, aes(x = estimator, y = estimate, fill = estimator)) +
  # FIX: outlier.shape = NA hides extreme heavy-tail outliers from the visual render.
  # This allows facet_wrap(scales = "free_y") to dynamically zoom in on the IQR/whiskers
  # without skewing the underlying statistical calculation of the boxplot itself.
  geom_boxplot(
    width = 0.68,
    outlier.shape = NA, 
    linewidth = 0.35
  ) +
  geom_hline(
    data = truth_df,
    aes(yintercept = truth_mu4),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.55,
    color = "firebrick"
  ) +
  geom_text(
    data = truth_df,
    aes(x = 1.5, y = truth_mu4, label = truth_lab),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 3.4,
    color = "firebrick"
  ) +
  facet_wrap(~ dist, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("Pairwise" = "#4C78A8", "Classical" = "#F58518")
  ) +
  labs(
    title = expression(paste("Pairwise vs Classical Estimation of ", mu[4])),
    subtitle = paste0("n = ", n, ", Monte Carlo replications = ", M, " (Visual outliers omitted for readability)"),
    x = NULL,
    y = expression(hat(mu)[4])
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  filename = "results/fig_pairwise_vs_classical.pdf",
  plot = p,
  width = 9.5,
  height = 7.0
)

ggsave(
  filename = "results/fig_pairwise_vs_classical.png",
  plot = p,
  width = 9.5,
  height = 7.0,
  dpi = 320
)

# -----------------
# Save outputs
# -----------------
write.csv(
  res,
  file = "results/fig_pairwise_vs_classical_raw.csv",
  row.names = FALSE
)

write.csv(
  summary_tab,
  file = "results/fig_pairwise_vs_classical_summary.csv",
  row.names = FALSE
)

print(summary_tab)