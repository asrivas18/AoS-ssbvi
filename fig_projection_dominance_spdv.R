# ==============================================================================
# fig_projection_dominance_spdv.R
# Projection Dominance Diagnostic (R_n) for SPDV_n (Figure 2)
# ==============================================================================
rm(list = ls())
source("functions_pairwise.R")
suppressPackageStartupMessages({ library(ggplot2); library(parallel); library(pbapply) })

set.seed(20260502)

n <- 100L
M <- 10000L
# Dense grid around the nu=4 threshold
nu_grid <- c(2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 16, 30)

cl <- makeCluster(max(1L, detectCores() - 1L))
clusterEvalQ(cl, { source("functions_pairwise.R") })

out_list <- lapply(nu_grid, function(nu) {
  distname <- paste0("t", nu)
  message(sprintf("Running nu = %.2f...", nu))
  clusterExport(cl, varlist = c("n", "distname"), envir = environment())
  
  results <- pblapply(seq_len(M), function(m) {
    x <- simulate_one_sample(n, distname, standardize_contam = FALSE)
    spdv_val <- spdv(x)
    # Projection variance for SPDV: tau_2^2 / n
    se_val <- se_spdv_projection(x, truncate_nonnegative = TRUE) 
    tau2_sq <- if (is.finite(se_val)) (se_val^2) * n else NA_real_
    list(spdv = spdv_val, tau2_sq = tau2_sq)
  }, cl = cl)
  
  spdv_vals <- sapply(results, `[[`, "spdv")
  tau2_sq_vals <- sapply(results, `[[`, "tau2_sq")
  
  valid_idx <- is.finite(spdv_vals) & is.finite(tau2_sq_vals)
  var_spdv <- stats::var(spdv_vals[valid_idx])
  avg_tau2_sq <- mean(tau2_sq_vals[valid_idx])
  
  R_n <- (avg_tau2_sq / n) / var_spdv
  data.frame(nu = nu, R_n = pmin(pmax(R_n, 0), 1.5))
})
stopCluster(cl)

res_rn <- do.call(rbind, out_list)

p <- ggplot(res_rn, aes(x = nu, y = R_n)) +
  geom_hline(yintercept = 1.0, color = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = 4.0, linetype = "dashed", color = "firebrick", linewidth = 0.85) +
  geom_line(color = "#4C78A8", linewidth = 1.3) +
  geom_point(fill = "#4C78A8", color = "white", shape = 21, size = 3.5, stroke = 1.1) +
  annotate("text", x = 4.2, y = 1.05, label = "nu = 4", color = "firebrick", angle = 90, vjust = -0.3) +
  scale_x_log10(breaks = c(3, 4, 5, 6, 8, 10, 16, 30), minor_breaks = NULL) +
  coord_cartesian(ylim = c(0, 1.25)) +
  labs(title = expression(paste("Projection Dominance Diagnostic for ", SPDV[n])),
       subtitle = "Linear Hoeffding projection vs. total variance (n = 100)",
       x = expression(paste("Degrees of freedom ", nu, " (log scale)")),
       y = expression(paste("Projection Ratio ", R[n]))) +
  theme_bw(base_size = 13)

ggsave("results/supp_projection_dominance_spdv.pdf", p, width=8.5, height=5.8)