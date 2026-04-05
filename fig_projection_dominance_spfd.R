# ==============================================================================
# fig_projection_dominance_spfd.R
# Publication-ready Figure: Projection Dominance Diagnostic (R_n) for SPFD_n
#
# R_n = (Projection variance) / (Total variance of SPFD_n)
# Shows breakdown of linear Hájek projection when 8th moment fails (nu <= 8)
#
# Matches Figure S3 / supp_projection_dominance_spfd.pdf in the supplement.
# ==============================================================================

rm(list = ls())
source("functions_pairwise.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(parallel)
  library(pbapply)      # for progress bar
})

set.seed(20260401)

# -----------------
# Parameters
# -----------------
n       <- 100L
M       <- 10000L     # High replication for stability under heavy tails
nu_grid <- c(4.5, 5, 6, 7, 7.5, 8, 8.5, 9, 10, 12, 16, 20, 30)

# -----------------
# Setup Cross-Platform Parallel Cluster
# -----------------
cores <- max(1L, detectCores() - 1L)
message(sprintf("Setting up parallel cluster with %d cores...", cores))

# Create cluster once to avoid massive overhead
cl <- makeCluster(cores)

# Force every worker core to load the pairwise functions
clusterEvalQ(cl, {
  source("functions_pairwise.R")
})

# -----------------
# Evaluator Function
# -----------------
evaluate_projection_dominance_spfd <- function(n, nu, M, cluster) {
  
  distname <- paste0("t", nu)
  message(sprintf("Running nu = %.2f (n = %d, M = %d)...", nu, n, M))
  
  # Export variables needed for this specific loop iteration
  clusterExport(cluster, varlist = c("n", "distname"), envir = environment())
  
  results <- pblapply(seq_len(M), function(m) {
    tryCatch({
      # Generate sample
      x <- simulate_one_sample(
        n = n,
        dist_name = distname,
        standardize_contam = FALSE
      )
      
      # Full statistic
      spfd_val <- spfd(x)
      
      # Projection-based standard error
      se_val <- se_spfd_projection(x, truncate_nonnegative = TRUE)
      tau4_sq <- if (is.finite(se_val)) (se_val^2) * n else NA_real_
      
      list(spfd = spfd_val, tau4_sq = tau4_sq, success = TRUE)
    }, error = function(e) {
      list(spfd = NA_real_, tau4_sq = NA_real_, success = FALSE)
    })
  }, cl = cluster)
  
  # Extract results
  spfd_vals    <- sapply(results, `[[`, "spfd")
  tau4_sq_vals <- sapply(results, `[[`, "tau4_sq")
  success_rate <- mean(sapply(results, `[[`, "success"))
  
  valid_idx <- !is.na(spfd_vals) & !is.na(tau4_sq_vals) & 
               is.finite(spfd_vals) & is.finite(tau4_sq_vals)
  
  var_spfd    <- if (sum(valid_idx) > 1) stats::var(spfd_vals[valid_idx]) else NA_real_
  avg_tau4_sq <- if (sum(valid_idx) > 0) mean(tau4_sq_vals[valid_idx]) else NA_real_
  
  R_n <- if (!is.na(var_spfd) && var_spfd > 0) {
    (avg_tau4_sq / n) / var_spfd
  } else NA_real_
  
  data.frame(
    nu           = nu,
    n            = n,
    M_total      = M,
    M_valid      = sum(valid_idx),
    success_rate = round(success_rate, 3),
    var_spfd     = var_spfd,
    avg_tau4_sq  = avg_tau4_sq,
    R_n          = R_n,
    stringsAsFactors = FALSE
  )
}

# -----------------
# Run simulation
# -----------------
out_list <- lapply(nu_grid, function(nu) {
  evaluate_projection_dominance_spfd(n = n, nu = nu, M = M, cluster = cl)
})

# Safely shut down the cluster
stopCluster(cl)

res_rn <- do.call(rbind, out_list)

# Cap for plotting (prevents chart from breaking on numeric instability below nu=6)
res_rn$R_n_plot <- pmin(pmax(res_rn$R_n, 0), 1.5)

# -----------------
# Create output directory
# -----------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# -----------------
# Publication-ready plot
# -----------------
p <- ggplot(res_rn, aes(x = nu, y = R_n_plot)) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 8.0, linetype = "dashed", color = "firebrick", linewidth = 0.85) +
  
  geom_line(color = "#4C78A8", linewidth = 1.3) +
  geom_point(fill = "#4C78A8", color = "white", shape = 21, size = 3.5, stroke = 1.1) +
  
  # Adjusted text annotations to fit the log-scaled x-axis
  annotate("text", x = 4.8, y = 1.12, 
           label = "Projection dominates\n(R_n ≈ 1)", 
           color = "#4C78A8", size = 4.2, fontface = "bold", hjust = 0) +
  
  annotate("text", x = 10, y = 0.25, 
           label = "Projection breaks down\n(8th moment fails)", 
           color = "firebrick", size = 4.2, hjust = 0, lineheight = 0.9) +
  
  annotate("text", x = 8.3, y = 1.05, 
           label = "nu = 8", color = "firebrick", size = 3.8, angle = 90, vjust = -0.3) +
  
  # Log scale perfectly stretches the critical nu in [4, 12] range
  scale_x_log10(
    breaks = c(4, 5, 6, 8, 10, 12, 16, 20, 30),
    minor_breaks = NULL
  ) +
  coord_cartesian(ylim = c(0, 1.25)) +
  
  labs(
    title = expression(paste("Projection Dominance Diagnostic for ", SPFD[n])),
    subtitle = bquote("Linear Hájek projection vs. total variance (n = 100, M = 10,000)"),
    x = expression(paste("Degrees of freedom ", nu, " (log scale)")),
    y = expression(paste("Projection Ratio ", R[n]))
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11.5, color = "grey30"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold", size = 14, margin = margin(r = 8)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.text = element_text(color = "black")
  )

# -----------------
# Save outputs
# -----------------
# Output names strictly match the LaTeX expectations
ggsave("results/supp_projection_dominance_spfd.pdf", plot = p, width = 8.5, height = 5.8, dpi = 300)
ggsave("results/supp_projection_dominance_spfd.png",  plot = p, width = 8.5, height = 5.8, dpi = 320)

write.csv(res_rn, "results/fig_projection_dominance_spfd_raw.csv", row.names = FALSE)

# Print summary
cat("\n=== SPFD Projection Dominance Diagnostic Complete ===\n")
print(res_rn[, c("nu", "M_valid", "success_rate", "R_n")], digits = 4)
cat("Files saved in 'results/' folder.\n")