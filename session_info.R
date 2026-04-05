# session_info.R
# ------------------------------------------------------------------------------
# Save R session information for reproducibility.
# ------------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE, recursive = TRUE)

outfile <- file.path("results", "sessionInfo.txt")

sink(outfile)
cat("Session information\n")
cat("===================\n\n")
cat("Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
sessionInfo()
sink()

message("Session information written to: ", outfile)