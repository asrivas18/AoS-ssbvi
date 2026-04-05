# ==============================================================================
# tab_mu4_comparison.R
# Generates Bias, Variance, and MSE for Table 3 in the supplement
# ==============================================================================
rm(list = ls())
source("functions_pairwise.R")
suppressPackageStartupMessages({ library(parallel); library(pbapply) })

set.seed(20260501)

# Design
n_grid <- c(50L, 100L)
dist_grid <- c("t6", "normal", "contam")
M <- 10000L

cores <- max(1L, detectCores() - 1L)
cl <- makeCluster(cores)
clusterEvalQ(cl, { source("functions_pairwise.R") })

# Run Simulation
results_list <- list()
for (dist_name in dist_grid) {
  for (n in n_grid) {
    message(sprintf("Running %s at n=%d...", dist_name, n))
    clusterExport(cl, varlist = c("n", "dist_name"), envir = environment())
    
    reps <- pblapply(seq_len(M), function(m) {
      x <- simulate_one_sample(n, dist_name, standardize_contam = TRUE)
      c(pair = mu4_pairwise(x), class = mu4_classical(x))
    }, cl = cl)
    
    mat <- do.call(rbind, reps)
    truth <- truth_mu4(dist_name, standardize_contam = TRUE)
    
    # Calculate Bias, Variance, MSE
    bias_pair <- mean(mat[, "pair"], na.rm=TRUE) - truth
    var_pair <- var(mat[, "pair"], na.rm=TRUE)
    mse_pair <- mean((mat[, "pair"] - truth)^2, na.rm=TRUE)
    
    bias_class <- mean(mat[, "class"], na.rm=TRUE) - truth
    var_class <- var(mat[, "class"], na.rm=TRUE)
    mse_class <- mean((mat[, "class"] - truth)^2, na.rm=TRUE)
    
    results_list[[paste(dist_name, n)]] <- data.frame(
      Distribution = ifelse(dist_name=="t6", "$t_6$", ifelse(dist_name=="normal", "Normal", "Contam. normal")),
      n = n,
      Estimator = c("Pairwise", "Classical"),
      Bias = c(bias_pair, bias_class),
      Variance = c(var_pair, var_class),
      MSE = c(mse_pair, mse_class)
    )
  }
}
stopCluster(cl)

# Format for LaTeX
final_tab <- do.call(rbind, results_list)
cat("\n\n--- LATEX TABLE ROWS ---\n")
for(i in 1:nrow(final_tab)) {
  cat(sprintf("%s & %d & %s & %.4f & %.4f & %.4f \\\\\n", 
              final_tab$Distribution[i], final_tab$n[i], final_tab$Estimator[i], 
              final_tab$Bias[i], final_tab$Variance[i], final_tab$MSE[i]))
}