# =============================================================================
# Comparison of glmnet (default), glmnet (proposed) and LARS on one train/test
# split. Used by the real data analysis (Section 4.3, Table 1).
#   run_lasso_comparison(): fits the three methods and returns test RMSE,
#                           computation time and SPE
#   count_sparsity():       number of nonzero coefficients / Jaccard indices
# =============================================================================
library(glmnet)
library(lars)
library(dplyr)
run_lasso_comparison <- function(y, x, ytest, xtest,
                                 thresh, nlambda,
                                 seed = 1,
                                 spe_points = 20,        # reference lambda points for SPE (k)
                                 spe_lambda_start = 0.001 # smallest reference lambda for SPE
) {
  
  # --- n and p are taken from the data ---
  n <- nrow(x)
  p <- ncol(x)
  
  has_warning_def  <- 0
  has_warning_lars <- 0
  has_warning_set  <- 0
  warning_msg_def  <- ""
  warning_msg_set  <- ""
  
  # =================================================================
  # 1. glmnet (default)
  # =================================================================
  if (!is.null(seed)) set.seed(seed)
  st <- Sys.time()
  lassocv <- cv.glmnet(x, y, family = "gaussian", alpha = 1)
  withCallingHandlers({
    model1 <- glmnet(x = x, y = y, family = "gaussian", alpha = 1)
  }, warning = function(w) {
    has_warning_def <<- 1
    warning_msg_def <<- paste(warning_msg_def, conditionMessage(w), sep = "; ")
    invokeRestart("muffleWarning")
  })
  gt <- Sys.time()
  
  lambda   <- lassocv$lambda.min
  pred1    <- predict(model1, newx = xtest, s = lambda)
  rmse_def <- c_rmse(pred1, ytest)
  time_def <- as.numeric(difftime(gt, st, units = "secs"))
  
  # =================================================================
  # 2. LARS
  #    mode = "step" if n == p, otherwise mode = "fraction"
  # =================================================================
  if (!is.null(seed)) set.seed(seed)
  lars_mode <- if (n == p) "fraction" else "step"
  if (n == p) {
    if (!is.null(seed)) set.seed(seed)
    st <- Sys.time()
    cvlars <- cv.lars(x, y, plot.it = FALSE, type = "lasso", mode = "fraction")
    withCallingHandlers({
      lmodel <- lars(x, y, type = "lasso", trace = FALSE)
    }, warning = function(w) {
      has_warning_lars <<- 1
      invokeRestart("muffleWarning")
    })
    gt <- Sys.time()
    beststep <- cvlars$index[which.min(cvlars$cv)]
    pred2    <- predict(lmodel, newx = xtest, mode = "fraction", s = beststep)
    rmsel    <- c_rmse(pred2$fit, ytest)
    timel    <- as.numeric(difftime(gt, st, units = "secs"))
    
  } else {
    if (!is.null(seed)) set.seed(seed)
    st <- Sys.time()
    cvlars <- cv.lars(x, y, plot.it = FALSE, type = "lasso", mode = "step", max.step = p,use.Gram=FALSE)
    withCallingHandlers({
      lmodel <- lars(x, y, type = "lasso", trace = FALSE)
    }, warning = function(w) {
      has_warning_lars <<- 1
      invokeRestart("muffleWarning")
    })
    gt <- Sys.time()
    beststep <- cvlars$index[which.min(cvlars$cv)]
    pred2    <- predict(lmodel, newx = xtest, mode = "step", s = beststep)
    rmsel    <- c_rmse(pred2$fit, ytest)
    timel    <- as.numeric(difftime(gt, st, units = "secs"))
  }
  
  # =================================================================
  # 3. Reference lambdas (candidate) for SPE and the exact LARS path
  # =================================================================
  lmax      <- max(model1$lambda)
  candidate <- exp(seq(from = log(spe_lambda_start), to = log(lmax),
                       length.out = spe_points))
  candidate[1] <- 0
  
  coef_lars <- predict.lars(lmodel, type = "coefficients",
                            mode = "lambda", s = candidate * sqrt(n))$coefficients
  coef_lars <- t(coef_lars)
  
  # --- SPE of glmnet (default) ---
  cg1          <- as.matrix(coef.glmnet(model1, s = candidate))[-1, ]
  rmsec1       <- apply(cg1 - coef_lars, 2, rmse)
  spe_def      <- mean(rmsec1)
  sd_spe_def   <- sd(rmsec1)
  
  # =================================================================
  # 4. glmnet (tuned): one configuration (thresh, nlambda)
  # =================================================================
  model1_lambda     <- model1$lambda
  len_model1_lambda <- length(model1_lambda)
  minn              <- min(model1_lambda)
  
  if (nlambda < len_model1_lambda) {
    # if nlambda is shorter than the default lambda sequence, use the default sequence
    lam <- model1_lambda
  } else {
    add <- seq(minn, 0, by = -(minn / (nlambda - len_model1_lambda)))[-1]
    lam <- c(model1_lambda, add)
  }
  
  if (!is.null(seed)) set.seed(seed)
  st <- Sys.time()
  cv_set <- cv.glmnet(x, y, family = "gaussian", alpha = 1,
                      thresh = thresh, lambda = lam, maxit = 1e+6)
  withCallingHandlers({
    glmnet_model <- glmnet(x, y, family = "gaussian", alpha = 1,
                           thresh = thresh, lambda = lam,maxit = 1e+6)
  }, warning = function(w) {
    has_warning_set <<- 1
    warning_msg_set <<- paste(warning_msg_set, conditionMessage(w), sep = "; ")
    invokeRestart("muffleWarning")
  })
  gt <- Sys.time()
  print(model1_lambda)
  print(lam)
  lambda_best <- cv_set$lambda.min
  time_glmnet <- as.numeric(difftime(gt, st, units = "secs"))
  
  # prediction error
  pred1       <- predict(glmnet_model, newx = xtest, s = lambda_best)
  rmse_glmnet <- c_rmse(pred1, ytest)
  
  # coefficient error (SPE)
  cg1            <- as.matrix(coef.glmnet(glmnet_model, s = candidate))[-1, ]
  rmsec1         <- apply(cg1 - coef_lars, 2, rmse)
  spe_glmnet     <- mean(rmsec1)
  sd_spe_glmnet  <- sd(rmsec1)
  
  # =================================================================
  # 5. Summary
  #    (LARS is the exact path, so its SPE is 0 by definition and not reported)
  # =================================================================
  summary_df <- data.frame(
    n = n, p = p,
    thresh = thresh, nlambda = nlambda,
    RMSE_pred_default = rmse_def,    time_default = time_def,
    RMSE_pred_lars    = rmsel,       time_lars    = timel,
    RMSE_pred_glmnet  = rmse_glmnet, time_glmnet  = time_glmnet,
    SPE_default = spe_def,    sd_SPE_default = sd_spe_def,
    SPE_glmnet  = spe_glmnet, sd_SPE_glmnet  = sd_spe_glmnet,
    best_lambda_default = lambda,
    best_lambda_glmnet  = lambda_best,
    best_step_lars = beststep,
    lars_mode = lars_mode,
    has_warning_def  = has_warning_def,
    has_warning_lars = has_warning_lars,
    has_warning_set  = has_warning_set
  )
  
  list(
    summary  = summary_df,
    models   = list(default = model1, lars = lmodel, tuned = glmnet_model),
    warnings = list(default = warning_msg_def, set = warning_msg_set)
  )
}


c_rmse<-function(pred,acc){
  l<-length(acc)
  rmse<-sqrt(sum((acc-pred)^2)/l)
  return(rmse)}


rmse<-function(x){
  l<-length(x)
  rmse<-sqrt(sum((x)^2)/l)
  return(rmse)
}

# ---- Sparse pattern counts (call after run_lasso_comparison) ----
count_sparsity <- function(res, tol = 0) {
  s   <- res$summary
  mod <- res$models
  
  # coefficients at each method's CV-selected point (intercept removed)
  b_def  <- as.vector(as.matrix(coef.glmnet(mod$default, s = s$best_lambda_default)))[-1]
  b_set  <- as.vector(as.matrix(coef.glmnet(mod$tuned,   s = s$best_lambda_glmnet)))[-1]
  b_lars <- as.vector(predict.lars(mod$lars, type = "coefficients",
                                   mode = "step", s = s$best_step_lars)$coefficients)
  
  nz_def  <- abs(b_def)  > tol
  nz_set  <- abs(b_set)  > tol
  nz_lars <- abs(b_lars) > tol
  p_all   <- length(b_def)
  
  data.frame(
    p_total      = p_all,
    nz_default   = sum(nz_def),
    nz_tuned     = sum(nz_set),
    nz_lars      = sum(nz_lars),
    zero_default = sum(!nz_def),
    zero_tuned   = sum(!nz_set),
    zero_lars    = sum(!nz_lars),
    # default vs tuned: how the supports differ
    common_def_set   = sum(nz_def & nz_set),
    only_default     = sum(nz_def & !nz_set),
    only_tuned       = sum(!nz_def & nz_set),
    jaccard_def_set  = sum(nz_def & nz_set) / max(1, sum(nz_def | nz_set)),
    # each glmnet vs LARS
    common_def_lars  = sum(nz_def & nz_lars),
    common_set_lars  = sum(nz_set & nz_lars),
    jaccard_def_lars = sum(nz_def & nz_lars) / max(1, sum(nz_def | nz_lars)),
    jaccard_set_lars = sum(nz_set & nz_lars) / max(1, sum(nz_set | nz_lars))
  )
}
