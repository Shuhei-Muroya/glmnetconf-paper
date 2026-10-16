# =============================================================================
# Appendix A: run glmnet over the 5 x 5 configuration grid for every dataset
# and record SPE (against LARS) and computation time (incl. CV).
# Output: A_configuration_effect/results/res_chunk_<j>.rds
#
# Run from the repository root.
# =============================================================================
library(doParallel)
library(foreach)
library(MASS)
library(glmnet)
library(Matrix)
library(lars)

source("A_configuration_effect/recipe.R")

# ---- settings ----
save_dir    <- file.path("A_configuration_effect", "results")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
start_chunk <- 1
core_n      <- 20
chunk_size  <- 20      # small chunks: 25 configurations per dataset
n_cand      <- 20      # number of reference lambdas for SPE
cand_min    <- 0.001

cl <- makeCluster(core_n)
registerDoParallel(cl)

master_recipe <- sensitivity_recipe_rep
n_total       <- nrow(master_recipe)
chunk_indices <- split(seq_len(n_total), ceiling(seq_len(n_total) / chunk_size))

clusterEvalQ(cl, {
  library(MASS); library(glmnet); library(Matrix); library(lars)
  source("R/data_generation.R")
})

clusterExport(cl, c("c_rmse", "rmse", "gen_data_flexible", "generate_sigma_wrapper",
                    "coef", "get_sparse_rates", "get_top_bottom_10",
                    "make_config_grid", "thresh_range", "n_grid_pts",
                    "nlam_min", "nlam_max_mult", "maxit_val",
                    "n_cand", "cand_min"))

cat("Starting from chunk", start_chunk, "of", length(chunk_indices), "\n")

for (j in start_chunk:length(chunk_indices)) {
  print(j)
  idx_j   <- chunk_indices[[j]]
  Plars_j <- master_recipe[idx_j, ]
  LN_j    <- nrow(Plars_j)
  
  result_df_j <- foreach(ln = 1:LN_j, .packages = c("MASS", "glmnet", "lars", "Matrix")) %dopar% {
    tryCatch({
      
      pattern <- Plars_j[ln, ]
      seed    <- as.numeric(pattern$seed)
      n       <- as.numeric(pattern$n)
      p       <- as.numeric(pattern$p)
      
      data_list_ln <- gen_data_flexible(pattern, seed = seed)
      x <- data_list_ln$x; y <- data_list_ln$y
      xtest <- data_list_ln$xtest; ytest <- data_list_ln$ytest
      
      # ---- default glmnet: base lambda sequence and the upper end (lmax) of the reference lambdas ----
      set.seed(seed)
      model1 <- glmnet(x = x, y = y, family = "gaussian", alpha = 1)
      model1_lambda     <- model1$lambda
      len_model1_lambda <- length(model1_lambda)
      minn <- min(model1_lambda)
      lmax <- max(model1_lambda)
      
      # ---- LARS (once per dataset) ----
      has_warning_lars <- 0
      set.seed(seed)
      lars_mode <- if (n != p) "fraction" else "step"
      st <- Sys.time()
      if (n != p) {
        cvlars <- cv.lars(x, y, plot.it = FALSE, type = "lasso", mode = lars_mode)
      } else {
        cvlars <- cv.lars(x, y, plot.it = FALSE, type = "lasso", mode = lars_mode, max.step = p)
      }
      withCallingHandlers({
        lmodel <- lars(x, y, type = "lasso", trace = FALSE)
      }, warning = function(w) {
        has_warning_lars <<- 1
        invokeRestart("muffleWarning")
      })
      timel    <- as.numeric(difftime(Sys.time(), st, units = "secs"))
      beststep <- cvlars$index[which.min(cvlars$cv)]
      pred2    <- predict(lmodel, newx = xtest, mode = lars_mode, s = beststep)
      rmsel    <- c_rmse(pred2$fit, ytest)
      coef_est_lars  <- as.vector(predict.lars(lmodel, type = "coefficients",
                                               mode = lars_mode, s = beststep)$coefficients)
      coef_rmse_lars <- rmse(coef_est_lars - data_list_ln$beta)
      
      candidate    <- exp(seq(from = log(cand_min), to = log(lmax), length.out = n_cand))
      candidate[1] <- 0
      coef_lars <- t(predict.lars(lmodel, type = "coefficients", mode = "lambda",
                                  s = candidate * sqrt(n))$coefficients)
      
      # ---- fit and evaluate glmnet once (used for the default call and for every grid point) ----
      run_glmnet <- function(extra_args) {
        has_warning_set <- 0
        set.seed(seed)
        st <- Sys.time()
        cv_set <- do.call(cv.glmnet, c(list(x = x, y = y, family = "gaussian", alpha = 1),
                                       extra_args))
        withCallingHandlers({
          fit <- do.call(glmnet, c(list(x = x, y = y, family = "gaussian", alpha = 1),
                                   extra_args))
        }, warning = function(w) {
          has_warning_set <<- 1
          invokeRestart("muffleWarning")
        })
        time_g      <- as.numeric(difftime(Sys.time(), st, units = "secs"))
        lambda_best <- cv_set$lambda.min
        
        pred_rmse <- c_rmse(predict(fit, newx = xtest, s = lambda_best), ytest)
        cg        <- as.matrix(coef.glmnet(fit, s = candidate))[-1, ]
        rmsec     <- apply(cg - coef_lars, 2, rmse)
        coef_est  <- as.vector(as.matrix(coef.glmnet(fit, s = lambda_best)))[-1]
        
        data.frame(
          RMSE_pred_glmnet = pred_rmse,
          time_glmnet      = time_g,
          SPE_glmnet       = mean(rmsec),
          sd_SPE_glmnet    = sd(rmsec),
          RMSE_coef_glmnet = rmse(coef_est - data_list_ln$beta),
          best_lambda      = lambda_best,
          has_warning_set  = has_warning_set
        )
      }
      
      # ---- default call (no arguments): check against the grid corner (1e-7, 100) ----
      res_default <- run_glmnet(list())
      res_default$config_id <- NA_integer_
      res_default$thresh    <- NA_real_
      res_default$nlambda   <- NA_real_
      res_default$is_default_call <- 1L
      
      # ---- run the whole grid ----
      config_grid <- make_config_grid(p, thresh_range, n_grid_pts, nlam_min, nlam_max_mult)
      
      res_list <- vector("list", nrow(config_grid))
      for (i in seq_len(nrow(config_grid))) {
        n_lambda <- config_grid$nlambda[i]
        thresh   <- config_grid$thresh[i]
        
        # lambda sequence: the default sequence extended by equally spaced values from its minimum to 0 (Sec. 3.2)
        if (n_lambda <= len_model1_lambda) {
          lam <- model1_lambda
        } else {
          add <- seq(minn, 0, by = -(minn / (n_lambda - len_model1_lambda)))[-1]
          lam <- c(model1_lambda, add)
        }
        
        res_i <- run_glmnet(list(thresh = thresh, lambda = lam, maxit = maxit_val))
        res_i$config_id <- config_grid$config_id[i]
        res_i$thresh    <- thresh
        res_i$nlambda   <- n_lambda
        res_i$is_default_call <- 0L
        res_list[[i]] <- res_i
      }
      
      glmnet_df <- do.call(rbind, c(list(res_default), res_list))
      glmnet_df$nlambda_over_p <- glmnet_df$nlambda / p
      
      lars_res <- list(
        n = n, p = p,
        case = pattern$case, pattern_id = pattern$pattern_id, rep_id = pattern$rep_id,
        sigma_struct = pattern$sigma_struct, rho = pattern$rho,
        beta_type = pattern$beta_type, beta_struct = pattern$beta_struct,
        x_dis_type = pattern$x_dis_type, snr = pattern$snr,
        error_sd = data_list_ln$noise_sd, gen_seed = seed,
        RMSE_pred_lars = rmsel, time_lars = timel,
        RMSE_coef_lars = coef_rmse_lars, has_warning_lars = has_warning_lars
      )
      
      list(lars_list = lars_res, glmnet_list = list(glmnet_df = glmnet_df))
      
    }, error = function(e) {
      list(error = TRUE, error_msg = conditionMessage(e),
           chunk = j, ln = ln, pattern = Plars_j[ln, ])
    })
  }
  
  saveRDS(result_df_j, file = file.path(save_dir, paste0("res_chunk_", j, ".rds")))
  
  error_tasks <- which(sapply(result_df_j, function(x) isTRUE(x$error)))
  if (length(error_tasks) > 0) {
    cat(sprintf("[%s] Chunk %d: %d task(s) failed (ln = %s)\n",
                Sys.time(), j, length(error_tasks), paste(error_tasks, collapse = ", ")))
  }
  cat(sprintf("[%s] Chunk %d/%d finished (%.1f%%)\n",
              Sys.time(), j, length(chunk_indices), 100 * j / length(chunk_indices)))
}

stopCluster(cl)