# =============================================================================
# Section 4.3.2 / Table 1: DeepCDR / GDSC, drug GDSC ID 1008
# (N = 475 cell lines, p = 697 genes, target = ln IC50)
# 10 random 9:1 train/test splits; glmnet (default), glmnet (proposed,
# T_hope = 10 s) and LARS. Output: 06_real_data/results/deepcdr_gdsc1008_Thope10.csv
#
# Run from the repository root (after prepare_deepcdr_gdsc.R).
# =============================================================================
source("R/lasso_comparison.R")   # run_lasso_comparison, c_rmse, rmse, count_sparsity
library(glmnetconf)
library(glmnet)
library(lars)

# ---- Settings (edit here only) ----
dataset_name <- "cancerdrug_deepcdr_gdsc_1008_N475_p697"
data_file    <- file.path("06_real_data", "data", "deepcdr_gdsc1008.csv")
train_ratio  <- 0.9
n_rep        <- 10
base_seed    <- 1
T_hope       <- 10
sparse_tol   <- 0                            # raise to 1e-8 if LARS looks noisy
out_dir      <- file.path("06_real_data", "results")

# ---- Data (target is the last column, predictors are all the others) ----
data        <- read.csv(file = data_file, check.names = FALSE)
target_col  <- tail(names(data), 1)          # taken from the data
first_x_col <- 1L
last_x_col  <- ncol(data) - 1L
n_all       <- nrow(data)

message("N = ", n_all, ", p = ", last_x_col, ", target = ", target_col)

# ---- Output paths ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv_name <- sprintf("deepcdr_gdsc1008_Thope%g.csv", T_hope)
csv_path <- file.path(out_dir, csv_name)

# ---- Repeat random splits ----
summary_all <- NULL
for (r in seq_len(n_rep)) {
  seed_r <- base_seed + r - 1
  set.seed(seed_r)
  
  train_idx  <- sample.int(n_all, floor(n_all * train_ratio))
  data_train <- data[train_idx, ]
  data_test  <- data[-train_idx, ]
  
  y_train <- data_train[[target_col]]
  x_train <- as.matrix(data_train[, first_x_col:last_x_col])
  y_test  <- data_test[[target_col]]
  x_test  <- as.matrix(data_test[, first_x_col:last_x_col])
  
  # configuration tuned on the training data only
  a <- tune_configuration(x_train, T_hope = T_hope)
  
  result <- run_lasso_comparison(y_train, x_train, y_test, x_test,
                                 nlambda = a$best_configuration$nlambda,
                                 thresh  = a$best_configuration$thresh,
                                 seed    = seed_r)
  
  sp <- count_sparsity(result, tol = sparse_tol)
  
  summary_all <- rbind(summary_all,
                       cbind(dataset = dataset_name, target = target_col,
                             rep = r, seed = seed_r, T_hope = T_hope,
                             result$summary, sp))
  
  # save after every repetition so that partial results survive a crash
  write.csv(summary_all, csv_path, row.names = FALSE)
  message("rep ", r, "/", n_rep, " done")
}

# ---- Quick look: mean and sd over repetitions ----
metric_cols <- c("RMSE_pred_default", "RMSE_pred_lars", "RMSE_pred_glmnet",
                 "time_default", "time_lars", "time_glmnet",
                 "SPE_default", "SPE_glmnet",
                 "nz_default", "nz_tuned", "nz_lars",
                 "jaccard_def_set", "jaccard_def_lars", "jaccard_set_lars")
print(rbind(mean = colMeans(summary_all[metric_cols]),
            sd   = apply(summary_all[metric_cols], 2, sd)))
