# =============================================================================
# Section 4.3: Table 1 (mean +- s.d. over the 10 random splits)
#
# Run from the repository root.
# =============================================================================
res_files <- c(`NRW.115`     = file.path("06_real_data", "results", "nrw115_clay_Thope5.csv"),
               `DeepCDR/GDSC` = file.path("06_real_data", "results", "deepcdr_gdsc1008_Thope10.csv"))

methods <- c(default = "glmnet (default)", glmnet = "glmnet (proposed)", lars = "LARS")

table1 <- do.call(rbind, lapply(names(res_files), function(ds) {
  d <- read.csv(res_files[[ds]])
  do.call(rbind, lapply(names(methods), function(m) {
    ms <- function(col, digits) {
      if (!col %in% names(d)) return("---")
      sprintf(paste0("%.", digits, "f +- %.", digits, "f"), mean(d[[col]]), sd(d[[col]]))
    }
    data.frame(Dataset = ds, Method = methods[[m]],
               RMSE = ms(paste0("RMSE_pred_", m), 3),
               SPE  = ms(paste0("SPE_", m), 4),
               Time = ms(paste0("time_", m), 1))
  }))
}))

print(table1, row.names = FALSE)
