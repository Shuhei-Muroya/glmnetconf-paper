# =============================================================================
# Section 4.1: glmnet (default) vs glmnet (proposed, T_hope = 20 s) vs LARS
#
# For each dataset: test RMSE, SPE against the LARS path, computation time
# (including cross-validation) and RMSE of the coefficients. The configuration
# of glmnet (proposed) is chosen by glmnetconf::tune_configuration().
# Output: 04_simulation/results/res_chunk_<j>.rds
#
# Run from the repository root.
# =============================================================================
library(doParallel)
library(foreach)
library(MASS)
library(glmnet)
library(Matrix)
library(lars)

source("04_simulation/recipe.R")
save_dir <- file.path("04_simulation", "results")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)


# tune_configuration() defaults: T_hope = 20 s, K = 1000 (not passed explicitly)
T_hope<-20
K<-1000 # number of sampled configurations

core_n <- 20
cl <- makeCluster(core_n)
registerDoParallel(cl)

# chunks
chunk_size <- 100
master_recipe <- simulation_recipe_rep
n_total <- nrow(master_recipe)
chunk_indices <- split(1:n_total, ceiling(seq_along(1:n_total) / chunk_size))

start_chunk <- 1

clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(Matrix)
  library(lars)
  library(glmnetconf)
  source("R/data_generation.R")
})

clusterExport(cl, c("c_rmse", "rmse", "gen_data_flexible", "generate_sigma_wrapper", 
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
      
      data_list_ln <- gen_data_flexible(pattern, seed = seed)
      # gen_data_flexible() also returns the true coefficients as data_list_ln$beta
      
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
      
      # --- error w.r.t. the true beta (glmnet default, coefficients at the CV-selected lambda) ---
      # note: coef() is masked by coef(p, sparse_ratio) defined in R/data_generation.R,
      #       so coef.glmnet() is called explicitly
      coef_est_def <- as.vector(as.matrix(coef.glmnet(model1, s = lambda)))[-1]
      coef_rmse_def <- rmse(coef_est_def - data_list_ln$beta)
      
      # LARS
      set.seed(seed)
      lars_mode <- if (n != p) "fraction" else "step"
      
      if (n != p) {
        st <- Sys.time()
        cvlars <- cv.lars(data_list_ln$x, data_list_ln$y, plot.it=FALSE, type="lasso", mode=lars_mode)
        withCallingHandlers({
          lmodel <- lars(data_list_ln$x, data_list_ln$y, type="lasso", trace=FALSE)
        }, warning = function(w) {
          has_warning_lars <<- 1
          invokeRestart("muffleWarning")
        })
        gt <- Sys.time()
        beststep <- cvlars$index[which.min(cvlars$cv)]
        pred2 <- predict(lmodel, newx=data_list_ln$xtest, mode=lars_mode, s=beststep)
        rmsel <- c_rmse(pred2$fit, data_list_ln$ytest)
        timel <- as.numeric(difftime(gt, st, units="secs"))
        
      } else {
        st <- Sys.time()
        cvlars <- cv.lars(data_list_ln$x, data_list_ln$y, plot.it=FALSE, type="lasso", mode=lars_mode, max.step=p)
        withCallingHandlers({
          lmodel <- lars(data_list_ln$x, data_list_ln$y, type="lasso", trace=FALSE)
        }, warning = function(w) {
          has_warning_lars <<- 1
          invokeRestart("muffleWarning")
        })
        gt <- Sys.time()
        beststep <- cvlars$index[which.min(cvlars$cv)]
        pred2 <- predict(lmodel, newx=data_list_ln$xtest, mode=lars_mode, s=beststep)
        rmsel <- c_rmse(pred2$fit, data_list_ln$ytest)
        timel <- as.numeric(difftime(gt, st, units="secs"))
      }
      
      # --- error w.r.t. the true beta (LARS, coefficients at the CV-selected step) ---
      coef_est_lars <- as.vector(predict.lars(lmodel, type="coefficients", mode=lars_mode, s=beststep)$coefficients)
      coef_rmse_lars <- rmse(coef_est_lars - data_list_ln$beta)
      
      # reference lambdas (candidate) for SPE (comparison with the LARS path)
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
        RMSE_coef_default=coef_rmse_def,   # error w.r.t. the true beta (glmnet default)
        RMSE_coef_lars=coef_rmse_lars,     # error w.r.t. the true beta (LARS)
        has_warning_def=has_warning_def,
        has_warning_lars=has_warning_lars
      )
      
      # configuration of glmnet (proposed) tuned by glmnetconf
      config <- tune_configuration(data_list_ln$x,seed=1,message = TRUE,line=FALSE)
      thresh_vec <- config$best_configuration$thresh
      plam <- floor(config$best_configuration$nlambda)
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
        
        # --- error w.r.t. the true beta (glmnet proposed, coefficients at lambda_best) ---
        # note: coef() is masked by coef(p, sparse_ratio) defined in R/data_generation.R,
        #       so coef.glmnet() is called explicitly
        coef_est_glmnet <- as.vector(as.matrix(coef.glmnet(glmnet_model, s = lambda_best)))[-1]
        coef_rmse_glmnet <- rmse(coef_est_glmnet - data_list_ln$beta)
        
        glmnet_df <- rbind(glmnet_df, c(n_lambda, thresh, pred_glmnet, time_glmnet, coef_accuracy_glmnet, coef_sd_glmnet, coef_rmse_glmnet, lambda_best, has_warning_set))
        warning_log[[i + 1]] <- list(
          type="set",
          n_lambda=Pattern_grid[i, 1],
          thresh=Pattern_grid[i, 2],
          msg=warning_msg_set
        )
      }
      
      glmnet_df <- as.data.frame(glmnet_df)
      colnames(glmnet_df) <- c("nlambda", "thresh", "RMSE_pred_glmnet", "time_glmnet",
                               "SPE_glmnet", "sd_SPE_glmnet", "RMSE_coef_glmnet",
                               "best_lambda", "has_warning_set")
      
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