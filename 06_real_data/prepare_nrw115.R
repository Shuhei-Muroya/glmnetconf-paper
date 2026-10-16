# =============================================================================
# Section 4.3.1: export the NRW.115 dataset from the LimeSoDa package
# (Schmidinger et al., 2025). N = 115 soil samples, p = 1686 mid-infrared
# bands; columns 1-3 are the targets (SOC, pH, Clay), predictors start at 4.
# Output: 06_real_data/data/nrw115.csv
#
# Run from the repository root.
# =============================================================================
library(LimeSoDa)

out_dir <- file.path("06_real_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(NRW.115$Dataset, file.path(out_dir, "nrw115.csv"), row.names = FALSE)
