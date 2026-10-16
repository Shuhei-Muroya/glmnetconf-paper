# =============================================================================
# Section 3.2.1: generation of the summary dataset
#
# For each recipe row: generate (X, y), fit glmnet (default) and LARS, then
# fit glmnet with the 9 configurations tau in {1e-7, 1e-8, 1e-9} x
# n_lambda in {floor(0.5p), p, floor(1.5p)} and record time and SPE against
# the LARS path. Results are saved chunk-wise (100 datasets per chunk) to
# 01_summary_dataset/results/<target>/res_chunk_<j>.rds.
#
# NOTE: this is computationally heavy (days on a many-core server). The
# outputs used in the paper are already stored in 01_summary_dataset/results/.
#
# Run from the repository root.
# =============================================================================
library(doParallel)
library(foreach)
library(MASS)
library(glmnet)
library(Matrix)
library(lars)

source("R/data_generation.R")
source("01_summary_dataset/recipe.R")

# Which part of the summary dataset to generate (one run per target)
target <- "high_snr"   # "xdist1_normal", "xdist2_lognormal", "xdist3_mixture",
                       # "xdist4_nonlinear_factor" or "high_snr"
# xdist1_normal was generated earlier (May 2026) with gen_data_flexible_v1();
# the other parts use gen_data_flexible().
recipe_list <- list(
  xdist1_normal           = master_recipe_xdis1,
  xdist2_lognormal        = master_recipe_xdis2,
  xdist3_mixture          = master_recipe_xdis3,
  xdist4_nonlinear_factor = master_recipe_xdis4,
  high_snr                = master_recipe_extra
)
gen_fun  <- if (target == "xdist1_normal") gen_data_flexible_v1 else gen_data_flexible
save_dir <- file.path("01_summary_dataset", "results", target)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# --- settings ---
core_n <- 20
cl <- makeCluster(core_n)
registerDoParallel(cl)

# chunks
chunk_size <- 100
master_recipe <- recipe_list[[target]]
n_total <- nrow(master_recipe)
chunk_indices <- split(1:n_total, ceiling(seq_along(1:n_total) / chunk_size))

start_chunk <- 1

clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(Matrix)
  library(lars)
  
  source("R/data_generation.R")
})

clusterExport(cl, c("c_rmse", "rmse", "gen_fun", "generate_sigma_wrapper", 
                    "coef", "get_sparse_rates", "get_top_bottom_10"))

cat("Resuming simulation from chunk", start_chunk, "of", length(chunk_indices), "...\n")

# --- main loop ---
for(j in start_chunk:length(chunk_indices)) {
  print(j)
  idx_j <- chunk_indices[[j]]
  Plars_j <- master_recipe[idx_j, ]
  LN_j <- nrow(Plars_j)
  
  result_df_j <- foreach(ln = 1:LN_j, .packages = c("MASS", "glmnet", "lars", "Matrix")) %dopar% {
    
    # record errors and continue
    tryCatch({
      
      pattern <- Plars_j[ln, ]
      seed <- as.numeric(pattern$seed)
      
      data_list_ln <- gen_fun(pattern, seed = seed)
      
      n <- as.numeric(pattern$n)
      p <- as.numeric(pattern$p)
      has_warning_def <- 0
      has_warning_lars <- 0
      warning_msg_def <- ""
      
      # glmnet (default)
      set.seed(seed)
      st <- Sys.time()
      lassocv <- cv.glmnet(data_list_ln$x, data_list_ln$y, family="gaussian", alpha=1)
      withCallingHandlers({
        model1 <- glmnet(x=data_list_ln$x, y=data_list_ln$y, family="gaussian", alpha=1)
      }, warning = function(w) {
        has_warning_def <<- 1
        warning_msg_def <<- paste(warning_msg_def, conditionMessage(w), sep="; ")
        invokeRestart("muffleWarning")
      })
      gt <- Sys.time()
      
      lambda <- lassocv$lambda.min
      pred1 <- predict(model1, newx=data_list_ln$xtest, s=lambda)
      rmse_def <- c_rmse(pred1, data_list_ln$ytest)
      time <- as.numeric(difftime(gt, st, units="secs"))
      
      # LARS
      set.seed(seed)
      if (n != p) {
        st <- Sys.time()
        cvlars <- cv.lars(data_list_ln$x, data_list_ln$y, plot.it=FALSE, type="lasso", mode="fraction")
        withCallingHandlers({
          lmodel <- lars(data_list_ln$x, data_list_ln$y, type="lasso", trace=FALSE)
        }, warning = function(w) {
          has_warning_lars <<- 1
          invokeRestart("muffleWarning")
        })
        gt <- Sys.time()
        beststep <- cvlars$index[which.min(cvlars$cv)]
        pred2 <- predict(lmodel, newx=data_list_ln$xtest, mode="fraction", s=beststep)
        rmsel <- c_rmse(pred2$fit, data_list_ln$ytest)
        timel <- as.numeric(difftime(gt, st, units="secs"))
        
      } else {
        st <- Sys.time()
        cvlars <- cv.lars(data_list_ln$x, data_list_ln$y, plot.it=FALSE, type="lasso", mode="step", max.step=p)
        withCallingHandlers({
          lmodel <- lars(data_list_ln$x, data_list_ln$y, type="lasso", trace=FALSE)
        }, warning = function(w) {
          has_warning_lars <<- 1
          invokeRestart("muffleWarning")
        })
        gt <- Sys.time()
        beststep <- cvlars$index[which.min(cvlars$cv)]
        pred2 <- predict(lmodel, newx=data_list_ln$xtest, mode="step", s=beststep)
        rmsel <- c_rmse(pred2$fit, data_list_ln$ytest)
        timel <- as.numeric(difftime(gt, st, units="secs"))
      }
      
      # reference lambdas (candidate) for SPE
      lmax <- max(model1$lambda)
      candidate <- exp(seq(from=log(0.001), to=log(lmax), length.out=20))
      candidate[1] <- 0
      
      cg1 <- as.matrix(coef.glmnet(model1, s=candidate))[-1, ]
      coef_lars <- predict.lars(lmodel, type="coefficients", mode="lambda", s=candidate*sqrt(n))$coefficients
      coef_lars <- t(coef_lars)
      
      rmsec1 <- apply(cg1 - coef_lars, 2, rmse)
      Krmsec1 <- mean(rmsec1)
      sdd <- sd(rmsec1)
      
      lars_res <- list(
        n=n, p=p,
        true_ev_top=data_list_ln$true_ev_top,
        true_ev_bottom=data_list_ln$true_ev_bottom,
        sample_ev_top=data_list_ln$sample_ev_top,
        sample_ev_bottom=data_list_ln$sample_ev_bottom,
        sigma_struct=pattern$sigma_struct,
        rho=pattern$rho,
        sprate=pattern$sprate,
        m=pattern$m,
        block_p=pattern$block_p,
        beta_type=pattern$beta_type,
        beta_struct=pattern$beta_struct,
        x_dis_type=pattern$x_dis_type,
        snr=pattern$snr,
        error_sd=data_list_ln$noise_sd,
        gen_seed=seed,
        RMSE_pred_default=rmse_def, time_default=time,
        RMSE_pred_lars=rmsel, time_lars=timel,
        SPE_default=Krmsec1, sd_SPE_default=sdd,
        has_warning_def=has_warning_def,
        has_warning_lars=has_warning_lars
      )
      
      thresh_vec <- c(1e-7, 1e-8, 1e-9)
      plam <- c(floor(p * 0.5), p, floor(p * 1.5))
      Pattern_grid <- expand.grid(plam=plam, thresh=thresh_vec)
      
      glmnet_df <- numeric()
      minn <- min(model1$lambda)
      model1_lambda <- model1$lambda
      len_model1_lambda <- length(model1$lambda)
      
      warning_log <- vector("list", nrow(Pattern_grid) + 1)
      warning_log[[1]] <- list(
        type="default",
        n_lambda=NA,
        thresh=NA,
        msg=warning_msg_def
      )
      
      for(i in 1:nrow(Pattern_grid)){
        has_warning_set <- 0
        warning_msg_set <- ""
        set.seed(seed)
        n_lambda <- Pattern_grid[i, 1]
        thresh <- Pattern_grid[i, 2]
        
        lam <- numeric()
        if(n_lambda < len_model1_lambda){
          lam <- model1_lambda
        } else {
          add <- seq(minn, 0, by=-(minn / (n_lambda - len_model1_lambda)))[-1]
          lam <- c(model1_lambda, add)
        }
        
        st <- Sys.time()
        cv_set <- cv.glmnet(data_list_ln$x, data_list_ln$y, family="gaussian", alpha=1, thresh=thresh, lambda=lam,maxit = 1e+07)
        withCallingHandlers({
          glmnet_model <- glmnet(data_list_ln$x, data_list_ln$y, family="gaussian", alpha=1, thresh=thresh, lambda=lam,maxit = 1e+07)
        }, warning = function(w) {
          has_warning_set <<- 1
          warning_msg_set <<- paste(warning_msg_set, conditionMessage(w), sep="; ")
          invokeRestart("muffleWarning")
        })
        gt <- Sys.time()
        lambda_best <- cv_set$lambda.min
        timek <- difftime(gt, st, units="secs")
        
        pred1 <- predict(glmnet_model, newx=data_list_ln$xtest, s=lambda_best)
        pred_glmnet <- c_rmse(pred1, data_list_ln$ytest)
        time_glmnet <- timek
        
        cg1 <- as.matrix(coef.glmnet(glmnet_model, s=candidate))[-1, ]
        rmsec1 <- apply(cg1 - coef_lars, 2, rmse)
        coef_accuracy_glmnet <- mean(rmsec1)
        coef_sd_glmnet <- sd(rmsec1)
        
        glmnet_df <- rbind(glmnet_df, c(n_lambda, thresh, pred_glmnet, time_glmnet, coef_accuracy_glmnet, coef_sd_glmnet, lambda_best, has_warning_set))
        warning_log[[i + 1]] <- list(
          type="set",
          n_lambda=Pattern_grid[i, 1],
          thresh=Pattern_grid[i, 2],
          msg=warning_msg_set
        )
      }
      
      glmnet_df <- as.data.frame(glmnet_df)
      colnames(glmnet_df) <- c("nlambda", "thresh", "RMSE_pred_glmnet", "time_glmnet", "SPE_glmnet", "sd_SPE_glmnet", "best_lambda", "has_warning_set")
      
      list(lars_list=lars_res, glmnet_list=list(glmnet_df=glmnet_df), warning_list=warning_log)
      
    }, error = function(e) {
      # record the failed task and skip it
      list(
        error = TRUE,
        error_msg = conditionMessage(e),
        chunk = j,
        ln = ln,
        pattern = Plars_j[ln, ]
      )
    })
  }
  
  # save each chunk
  saveRDS(result_df_j, file=file.path(save_dir, paste0("res_chunk_", j, ".rds")))
  
  # report failed tasks
  error_tasks <- which(sapply(result_df_j, function(x) isTRUE(x$error)))
  if(length(error_tasks) > 0) {
    cat(sprintf("[%s] Chunk %d: %d task(s) failed (ln = %s)\n",
                Sys.time(), j, length(error_tasks), paste(error_tasks, collapse=", ")))
  }
  
  cat(sprintf("[%s] Chunk %d/%d finished (%.1f%%)\n",
              Sys.time(), j, length(chunk_indices), 100 * j / length(chunk_indices)))
}

stopCluster(cl)