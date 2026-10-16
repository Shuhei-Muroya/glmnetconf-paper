# =============================================================================
# Section 3.2.1: assemble the summary dataset used to train glmnet-MLP
#
# Reads the chunked outputs of run_generation.R and writes one row per
# (dataset, configuration):
#   (N, p, gamma_1..gamma_10, gamma_-1..gamma_-10, n_lambda, tau, T_glmnet, SPE)
# = 63,000 datasets x 9 configurations = 567,000 rows.
# The eigenvalues are those of the sample covariance matrix (sample_ev_*).
#
# Run from the repository root.
# =============================================================================
library(dplyr)
library(purrr)
source("R/load_results.R")

# ---- Settings ----
result_dirs <- file.path("01_summary_dataset", "results",
                         c("xdist1_normal", "xdist2_lognormal", "xdist3_mixture",
                           "xdist4_nonlinear_factor", "high_snr"))
out_file    <- file.path("01_summary_dataset", "summary_dataset.csv.gz")

x_cols <- c("n", "p",
            paste0("sample_ev_top", 1:10),
            paste0("sample_ev_bottom", 1:10),
            "nlambda", "thresh")
y_cols <- c("time_glmnet", "SPE_glmnet")

# ---- Combine ----
summary_dataset <- map_dfr(result_dirs, function(d) {
  process_simulation_results(d)$all_combined
}) %>%
  select(all_of(c(x_cols, y_cols))) %>%
  mutate(across(everything(), as.numeric))

cat("rows:", nrow(summary_dataset), "\n")   # expected: 567000

write.csv(summary_dataset, gzfile(out_file), row.names = FALSE)
