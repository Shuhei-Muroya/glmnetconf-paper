# =============================================================================
# Data generation for all simulation studies (Section 3.2.1, Appendix B.1)
#
# gen_data_flexible(recipe, seed) generates (X, y) for one row of a recipe:
#   x_dis_type 1: multivariate normal     2: log-normal
#              3: mixture of normals      4: nonlinear factor model
#              5: multivariate t (df = 5; used only in Section 4.1 and App. D)
#   beta_struct 1: nonzero coefficients = 1
#               2: nonzero coefficients = sign * (1 + |z|)
#   beta_type   1..5: number of nonzeros = log(p), sqrt(p), 0.1p, 0.3p, 0.5p
# It also returns the top/bottom 10 eigenvalues of the true and sample
# covariance matrices, which are the input features of glmnet-MLP.
#
# gen_data_flexible_v1() is the earlier revision used for xdist1_normal.
#
# All scripts are run with the repository root as the working directory.
# =============================================================================
source("R/covariance.R")

library(MASS)
library(dplyr)

generate_sigma_wrapper <- function(p, recipe, seed) {
  struct <- recipe$sigma_struct
  
  # values are taken from one row of the recipe;
  # columns that are NA are used only where needed
  if (struct == "compound") {
    return(generate_compound_sigma(p, rho = recipe$rho))
  } else if (struct == "AR") {
    return(generate_ar_sigma(p, rho = recipe$rho))
  } else if (struct == "random_sp") {
    return(random_matrix_sigma_sparse(p, sprate = recipe$sprate, seed = seed))
  } else if (struct == "ba_hub") {
    return(generate_ba_hub_sigma(p, m = recipe$m, seed = seed))
  } else if (struct == "block") {
    # block_p is the number of blocks
    return(generate_block_sigma(p, n_blocks = recipe$block_p, seed = seed))
  }
}




get_top_bottom_10 <- function(ev) {
  list(
    top10 = ev[1:min(10, length(ev))],
    bottom10 = rev(tail(ev, min(10, length(ev)))) # smallest first
  )
}



## Coefficient vectors
coef1 <- function(p, sparse_ratio = 0.9) {
  beta <- rep(0, p)
  n_active <- round(p * (1 - sparse_ratio))
  
  if (n_active > 0) {
    active_idx <- sample(1:p, n_active)
    beta[active_idx] <- 1
  }
  return(beta)
}

coef <- function(p, sparse_ratio = 0.9) {
  beta <- rep(0, p)
  n_active <- round(p * (1 - sparse_ratio))
  
  if (n_active > 0) {
    active_idx <- sample(1:p, n_active)
    
    # 4. Sign * (1 + |z|)
    S <- sample(c(-1, 1), n_active, replace = TRUE)
    z <- rnorm(n_active)
    beta[active_idx] <- S * (1 + abs(z))
  }
  return(beta)
}

# Five sparsity levels of beta (fraction of zeros) for a given p
get_sparse_rates <- function(p) {
  # 1. log(p) nonzero coefficients
  s_log <- max(1, round(log(p))) 
  rate_log <- 1 - (s_log / p)
  
  # 2. sqrt(p) nonzero coefficients
  s_sqrt <- max(1, round(sqrt(p)))
  rate_sqrt <- 1 - (s_sqrt / p)
  
  # 3-5. 10%, 30% and 50% nonzero coefficients
  rate_10 <- 0.90
  rate_30<-0.7
  rate_50<-0.5
  
  # keep all rates below 1 (for very small p)
  return(pmin(0.99, c(rate_log, rate_sqrt, rate_10,rate_30,rate_50)))
}


gen_data_flexible <- function(recipe, seed, n_test = 100) {
  set.seed(seed)
  
  # 1. Parameters
  n <- as.numeric(recipe$n)
  p <- as.numeric(recipe$p)
  snr <- as.numeric(recipe$snr)
  x_dis <- as.numeric(recipe$x_dis_type)
  beta_struct <- as.numeric(recipe$beta_struct)
  
  # 2. Covariance matrix Sigma (of the latent Gaussian Z; used for generation and true_ev)
  sigma <- generate_sigma_wrapper(p, recipe, seed)
  
  # 3. Coefficient vector
  sparse_rates <- get_sparse_rates(p)
  target_rate  <- sparse_rates[as.numeric(recipe$beta_type)]
  if (beta_struct == 1) {
    beta <- coef1(p, sparse_ratio = target_rate)
  } else {
    beta <- coef(p, sparse_ratio = target_rate)
  }
  
  # 4. Design matrix X
  N_total <- n + n_test
  mu <- rep(0, p)
  
  if (x_dis == 1) {
    # multivariate normal
    X_all <- mvrnorm(N_total, mu, sigma)
  } else if (x_dis == 2) { # lognormal
    Z     <- mvrnorm(N_total, mu, sigma)
    X_all <- exp(Z)
    # single component, so no empty-component issue;
    # column scaling is done below for all designs
    
  } else if (x_dis == 3) { # mixture of normals (bimodal)
    d_mix <- 0.5
    
    # mean of the second component: shifted by d_mix along the first principal direction of Sigma
    ev      <- eigen(sigma)
    v1      <- ev$vectors[, 1]
    lambda1 <- ev$values[1]
    mu2     <- d_mix * sqrt(lambda1) * v1
    
    # equal-weight Bernoulli assignment (eps = 0.5)
    eps  <- 0.5
    comp <- rbinom(N_total, 1, eps)
    
    # safeguard: at least one observation in each component
    if (sum(comp == 0) == 0) comp[1]        <- 0L
    if (sum(comp == 1) == 0) comp[N_total]  <- 1L
    
    n_comp1 <- sum(comp == 0)
    n_comp2 <- sum(comp == 1)
    
    X_comp1 <- mvrnorm(n_comp1, mu = rep(0, p), Sigma = sigma)
    X_comp2 <- mvrnorm(n_comp2, mu = mu2,       Sigma = sigma)
    
    X_all              <- matrix(NA, N_total, p)
    X_all[comp == 0, ] <- X_comp1
    X_all[comp == 1, ] <- X_comp2
    
  } else if (x_dis == 4) { # nonlinear factor model (shared factors)
    gamma <- 0.5  
    m_nl  <- max(2L, floor(p / 5))
    
    # --- linear part ---
    X_linear <- mvrnorm(N_total, mu, sigma)
    
    # --- factor assignment ---
    factor_assign <- sample(seq_len(m_nl), p, replace = TRUE)
    
    # --- factor matrix: N_total x m_nl ---
    C_matrix <- matrix(rnorm(N_total * m_nl), N_total, m_nl)
    
    # --- nonlinear function chosen at random for each variable j ---
    h_choice <- sample(c("sq", "tanh", "exp"), p, replace = TRUE)
    
    # --- nonlinear part (vectorised) ---
    X_nonlinear <- matrix(0, N_total, p)
    
    idx_sq   <- which(h_choice == "sq")
    idx_tanh <- which(h_choice == "tanh")
    idx_exp  <- which(h_choice == "exp")
    
    if (length(idx_sq)   > 0)
      X_nonlinear[, idx_sq]   <- C_matrix[, factor_assign[idx_sq],   drop = FALSE]^2
    if (length(idx_tanh) > 0)
      X_nonlinear[, idx_tanh] <- tanh(C_matrix[, factor_assign[idx_tanh], drop = FALSE])
    if (length(idx_exp)  > 0)
      X_nonlinear[, idx_exp]  <- exp(C_matrix[,  factor_assign[idx_exp],  drop = FALSE])
    
    X_all <- X_linear + gamma * X_nonlinear
    # column scaling is done below (sweep)
  } else {
    nu <- 5  # degrees of freedom
    
    Z <- mvrnorm(N_total, mu, sigma)
    W <- rchisq(N_total, df = nu)
    X_all <- Z / sqrt(W / nu)
  }
  
  # --- all designs: scale columns to unit variance (no centring) ---
  # so that SNR, eigenvalue features and SPE are comparable across designs;
  # centring is done inside glmnet.
  col_sd <- apply(X_all, 2, sd)
  col_sd[col_sd == 0] <- 1                    # guard against constant columns (division by zero)
  X_all  <- sweep(X_all, 2, col_sd, FUN = "/")
  
  # training / test split
  x      <- X_all[1:n, ]
  xtest  <- X_all[(n + 1):N_total, ]
  
  # 5. Response y with noise determined by the SNR
  # SNR = Var(X beta) / sigma^2, where Var(X beta) = beta' S beta with S the sample
  # covariance of the scaled design (beta' Sigma beta for Gaussian designs;
  # the actual covariance after the transformation for non-Gaussian designs).
  signal     <- x %*% beta
  signal_var <- as.numeric(var(as.vector(signal)))
  noise_sd   <- sqrt(signal_var / snr)
  errors     <- rnorm(N_total, 0, noise_sd)
  
  y      <- x     %*% beta + errors[1:n]
  ytest  <- xtest %*% beta + errors[(n + 1):N_total]
  
  # 6. Eigenvalue features
  # eigenvalues of the true covariance matrix (Sigma, latent space)
  true_ev_all <- eigen(sigma, only.values = TRUE)$values
  true_ev     <- get_top_bottom_10(true_ev_all)
  
  # eigenvalues of the sample covariance matrix (of the scaled X)
  sample_sigma  <- (t(x) %*% x) / n
  sample_ev_all <- eigen(sample_sigma, only.values = TRUE)$values
  sample_ev     <- get_top_bottom_10(sample_ev_all)
  
  # return everything as a list
  return(list(
    x = x, y = y, xtest = xtest, ytest = ytest,
    beta = beta,
    true_ev_top = true_ev$top10,
    true_ev_bottom = true_ev$bottom10,
    sample_ev_top = sample_ev$top10,
    sample_ev_bottom = sample_ev$bottom10,
    snr = snr,
    noise_sd = noise_sd
  ))
}


# -----------------------------------------------------------------------------
# gen_data_flexible_v1(): earlier revision of gen_data_flexible(), used only for
# the multivariate-normal part of the summary dataset
# (01_summary_dataset/results/xdist1_normal, generated in May 2026).
# Differences for x_dis_type = 1: no column scaling of X, and the noise
# variance is beta' Sigma beta / SNR. (Its other x_dis_type branches were not
# used for any result.)
# -----------------------------------------------------------------------------

gen_data_flexible_v1 <- function(recipe, seed, n_test = 100) {
  set.seed(seed)
  
  # 1. Parameters
  n <- as.numeric(recipe$n)
  p <- as.numeric(recipe$p)
  snr <- as.numeric(recipe$snr)
  x_dis <- as.numeric(recipe$x_dis_type)
  beta_struct<-as.numeric(recipe$beta_struct)
  
  # 2. Covariance matrix Sigma
  sigma <- generate_sigma_wrapper(p, recipe, seed)
  
  # 3. Coefficient vector
  # beta_type 1-5: log(p), sqrt(p), 10%, 30%, 50% nonzero coefficients
  sparse_rates <- get_sparse_rates(p)
  target_rate <- sparse_rates[as.numeric(recipe$beta_type)]
  if(beta_struct==1){
    beta <- coef1(p, sparse_ratio = target_rate)
  }else{
    beta <- coef(p, sparse_ratio = target_rate)
  }
  
  # 4. Design matrix X
  N_total <- n + n_test
  mu <- rep(0, p)
  
  if (x_dis == 1) {
    # multivariate normal
    X_all <- mvrnorm(N_total, mu, sigma)
  }else if (x_dis == 2) {#lognormal
    Z <- mvrnorm(N_total, mu, sigma)
    X_all <- exp(Z)
  } else if (x_dis == 3) {#WeightedLAD-LASSO
    # epsilon-contamination: (1-eps)*N(mu, Sigma) + eps*N(mu, k^2 * Sigma)
    eps <- 0.05
    k   <- 5
    
    n_contam <- rbinom(1, N_total, eps)
    n_clean  <- N_total - n_contam
    
    # keep the order of random number calls identical to the original code
    X_clean <- mvrnorm(n_clean, mu, sigma)
    if (n_contam > 0) {
      X_contam <- mvrnorm(n_contam, mu, k^2 * sigma)
    }
    row_idx <- sample(N_total)
    
    X_all <- matrix(NA, N_total, p)
    X_all[row_idx[1:n_clean], ] <- X_clean
    if (n_contam > 0) {
      X_all[row_idx[(n_clean + 1):N_total], ] <- X_contam
    }
  }else {
    # multivariate t (3 degrees of freedom):
    # normal draws divided by sqrt(chi-sq / df)
    Z <- mvrnorm(N_total, mu, sigma)
    W <- rchisq(N_total, df = 3)
    X_all <- Z / sqrt(W / 3)
  }
  
  # training / test split
  x      <- X_all[1:n, ]
  xtest  <- X_all[(n + 1):N_total, ]
  
  # 5. Response y with noise determined by the SNR (beta' Sigma beta / SNR)
  signal_var <- as.numeric(t(beta) %*% sigma %*% beta)
  noise_sd   <- sqrt(signal_var / snr)
  errors     <- rnorm(N_total, 0, noise_sd)
  
  y      <- x %*% beta + errors[1:n]
  ytest  <- xtest %*% beta + errors[(n + 1):N_total]
  
  # 6. Eigenvalue features
  # eigenvalues of the true covariance matrix
  true_ev_all <- eigen(sigma, only.values = TRUE)$values
  true_ev     <- get_top_bottom_10(true_ev_all)
  
  # eigenvalues of the sample covariance matrix
  sample_sigma  <- (t(x) %*% x) / n
  sample_ev_all <- eigen(sample_sigma, only.values = TRUE)$values
  sample_ev     <- get_top_bottom_10(sample_ev_all)
  
  # return everything as a list
  return(list(
    x = x, y = y, xtest = xtest, ytest = ytest,
    beta = beta,
    true_ev_top = true_ev$top10,
    true_ev_bottom = true_ev$bottom10,
    sample_ev_top = sample_ev$top10,
    sample_ev_bottom = sample_ev$bottom10,
    snr = snr,
    noise_sd=noise_sd
  ))
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