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
  
  # recipe (1行分のデータフレーム) から各値を取り出す
  # NAが入っている可能性があるので、必要な時だけ使う
  if (struct == "compound") {
    return(generate_compound_sigma(p, rho = recipe$rho))
  } else if (struct == "AR") {
    return(generate_ar_sigma(p, rho = recipe$rho))
  } else if (struct == "random_sp") {
    return(random_matrix_sigma_sparse(p, sprate = recipe$sprate, seed = seed))
  } else if (struct == "ba_hub") {
    return(generate_ba_hub_sigma(p, m = recipe$m, seed = seed))
  } else if (struct == "block") {
    # block_p はここでは「ブロック数」として扱う
    return(generate_block_sigma(p, n_blocks = recipe$block_p, seed = seed))
  }
}

# random
#sigma_list<-list(f1=sigma1,f2=sigma2,f3=random_matrix_sigma_sparse,f4=generate_ba_hub_sigma,f5=generate_block_sigma)



get_top_bottom_10 <- function(ev) {
  list(
    top10 = ev[1:min(10, length(ev))],
    bottom10 = rev(tail(ev, min(10, length(ev)))) # 小さい順
  )
}



##Coef 係数の関数
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

# p に基づいて 4 種類の sparse_ratio を返す関数
get_sparse_rates <- function(p) {
  # 例: log(p) / p non 0
  s_log <- max(1, round(log(p))) 
  rate_log <- 1 - (s_log / p)
  
  # 2. 中間設定1 (O(sqrt(p)))
  s_sqrt <- max(1, round(sqrt(p)))
  rate_sqrt <- 1 - (s_sqrt / p)
  
  # 3. 中間設定2 (1/20 = 5% が非ゼロ)
  rate_10 <- 0.90
  rate_30<-0.7
  rate_50<-0.5
  
  # 念のため、すべて 1 未満になるように調整（pが極端に小さい場合用）
  return(pmin(0.99, c(rate_log, rate_sqrt, rate_10,rate_30,rate_50)))
}


gen_data_flexible <- function(recipe, seed, n_test = 100) {
  set.seed(seed)
  
  # 1. パラメータの抽出
  n <- as.numeric(recipe$n)
  p <- as.numeric(recipe$p)
  snr <- as.numeric(recipe$snr)
  x_dis <- as.numeric(recipe$x_dis_type)
  beta_struct <- as.numeric(recipe$beta_struct)
  
  # 2. 共分散行列 Sigma の生成（潜在ガウス空間 Z の共分散：生成と true_ev に使用）
  sigma <- generate_sigma_wrapper(p, recipe, seed)
  
  # 3. 係数ベクトルの生成
  sparse_rates <- get_sparse_rates(p)
  target_rate  <- sparse_rates[as.numeric(recipe$beta_type)]
  if (beta_struct == 1) {
    beta <- coef1(p, sparse_ratio = target_rate)
  } else {
    beta <- coef(p, sparse_ratio = target_rate)
  }
  
  # 4. デザイン行列 X の生成
  N_total <- n + n_test
  mu <- rep(0, p)
  
  if (x_dis == 1) {
    # 正規分布
    X_all <- mvrnorm(N_total, mu, sigma)
  } else if (x_dis == 2) { # lognormal
    Z     <- mvrnorm(N_total, mu, sigma)
    X_all <- exp(Z)
    # n_comp=0 の心配なし（全サンプル同一成分）
    # 列標準化は共通処理で実施
    
  } else if (x_dis == 3) { # 混合正規（bimodal）
    d_mix <- 0.5
    
    # 第2成分の平均: Σの第1主成分方向にマハラノビス距離 d_mix
    ev      <- eigen(sigma)
    v1      <- ev$vectors[, 1]
    lambda1 <- ev$values[1]
    mu2     <- d_mix * sqrt(lambda1) * v1
    
    # 等重み（eps=0.5）でベルヌーイ割り当て
    eps  <- 0.5
    comp <- rbinom(N_total, 1, eps)
    
    # ② 保険: n_comp1=0 または n_comp2=0 のとき最低1件ずつ保証
    if (sum(comp == 0) == 0) comp[1]        <- 0L
    if (sum(comp == 1) == 0) comp[N_total]  <- 1L
    
    n_comp1 <- sum(comp == 0)
    n_comp2 <- sum(comp == 1)
    
    X_comp1 <- mvrnorm(n_comp1, mu = rep(0, p), Sigma = sigma)
    X_comp2 <- mvrnorm(n_comp2, mu = mu2,       Sigma = sigma)
    
    X_all              <- matrix(NA, N_total, p)
    X_all[comp == 0, ] <- X_comp1
    X_all[comp == 1, ] <- X_comp2
    
  } else if (x_dis == 4) { # 非線形因子モデル（因子共有型）
    gamma <- 0.5  
    m_nl  <- max(2L, floor(p / 5))
    
    # --- 線形部分 ---
    X_linear <- mvrnorm(N_total, mu, sigma)
    
    # --- 因子割り当て ---
    factor_assign <- sample(seq_len(m_nl), p, replace = TRUE)
    
    # --- 因子行列: N_total × m_nl ---
    C_matrix <- matrix(rnorm(N_total * m_nl), N_total, m_nl)
    
    # --- 非線形変換関数をランダム選択（変数 j ごと）---
    h_choice <- sample(c("sq", "tanh", "exp"), p, replace = TRUE)
    
    # --- ③ 非線形部分（ベクトル化・for ループなし）---
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
    # ↑ 列標準化は共通処理（sweep）で実施
  } else {
    nu <- 5  # degrees of freedom
    
    Z <- mvrnorm(N_total, mu, sigma)
    W <- rchisq(N_total, df = nu)
    X_all <- Z / sqrt(W / nu)
  }
  
  # --- 全パターン共通：列スケール標準化（分散1へ。中心化はしない）---
  # 理由: SNR・固有値特徴量・SPE をパターン間で揃えるため。
  #       中心化は glmnet が内部で行うため生成段階では行わない。
  col_sd <- apply(X_all, 2, sd)
  col_sd[col_sd == 0] <- 1                    # 退化列の保険（ゼロ割回避）
  X_all  <- sweep(X_all, 2, col_sd, FUN = "/")
  
  # 訓練データとテストデータの分割
  x      <- X_all[1:n, ]
  xtest  <- X_all[(n + 1):N_total, ]
  
  # 5. 応答変数 y の生成（SNRに基づくノイズ付加）
  # SNR = Var(Xβ) / σ^2 で定義。
  # Var(Xβ) は標準化後デザインの標本共分散 S を用いた β'Sβ に等しい
  # （ガウス設計では β'Σβ に一致、非ガウス設計では変換後の実共分散を反映）。
  signal     <- x %*% beta
  signal_var <- as.numeric(var(as.vector(signal)))
  noise_sd   <- sqrt(signal_var / snr)
  errors     <- rnorm(N_total, 0, noise_sd)
  
  y      <- x     %*% beta + errors[1:n]
  ytest  <- xtest %*% beta + errors[(n + 1):N_total]
  
  # 6. 特徴量（固有値）の抽出
  # 真の共分散行列の固有値（Σ 由来。生成と同じ潜在空間）
  true_ev_all <- eigen(sigma, only.values = TRUE)$values
  true_ev     <- get_top_bottom_10(true_ev_all)
  
  # サンプル共分散行列の固有値（標準化後の X 由来）
  sample_sigma  <- (t(x) %*% x) / n
  sample_ev_all <- eigen(sample_sigma, only.values = TRUE)$values
  sample_ev     <- get_top_bottom_10(sample_ev_all)
  
  # 全てをリストにまとめて返す
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
  
  # 1. パラメータの抽出
  n <- as.numeric(recipe$n)
  p <- as.numeric(recipe$p)
  snr <- as.numeric(recipe$snr)
  x_dis <- as.numeric(recipe$x_dis_type)
  beta_struct<-as.numeric(recipe$beta_struct)
  
  # 2. 共分散行列 Sigma の生成
  # 以前作成した generate_sigma_wrapper を利用
  sigma <- generate_sigma_wrapper(p, recipe, seed)
  
  # 3. 係数ベクトルの生成
  # beta_type (1:log, 2:sqrt, 3:5%, 4:10%, 5:20%,6:30%,7:50%) に応じて比率を決定
  sparse_rates <- get_sparse_rates(p)
  target_rate <- sparse_rates[as.numeric(recipe$beta_type)]
  if(beta_struct==1){
    beta <- coef1(p, sparse_ratio = target_rate)
  }else{
    beta <- coef(p, sparse_ratio = target_rate)
  }
  
  # 4. デザイン行列 X の生成
  N_total <- n + n_test
  mu <- rep(0, p)
  
  if (x_dis == 1) {
    # 正規分布
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
    
    # 乱数の呼び出し順序を元コードと完全に一致させる
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
    # t分布 (自由度3) などのバリエーション
    # 正規乱数を sqrt(chi-sq/df) で割る生成法
    Z <- mvrnorm(N_total, mu, sigma)
    W <- rchisq(N_total, df = 3)
    X_all <- Z / sqrt(W / 3)
  }
  
  # 訓練データとテストデータの分割
  x      <- X_all[1:n, ]
  xtest  <- X_all[(n + 1):N_total, ]
  
  # 5. 応答変数 y の生成 (SNRに基づくノイズ付加)
  signal_var <- as.numeric(t(beta) %*% sigma %*% beta)
  noise_sd   <- sqrt(signal_var / snr)
  errors     <- rnorm(N_total, 0, noise_sd)
  
  y      <- x %*% beta + errors[1:n]
  ytest  <- xtest %*% beta + errors[(n + 1):N_total]
  
  # 6. 特徴量（固有値）の抽出
  # 真の共分散行列の固有値
  true_ev_all <- eigen(sigma, only.values = TRUE)$values
  true_ev     <- get_top_bottom_10(true_ev_all)
  
  # サンプル共分散行列の固有値 (P > N の場合に注意)
  # メモリ節約のため sample_sigma を直接作らず tcrossprod 等を使うのも手
  sample_sigma  <- (t(x) %*% x) / n
  sample_ev_all <- eigen(sample_sigma, only.values = TRUE)$values
  sample_ev     <- get_top_bottom_10(sample_ev_all)
  
  # 全てをリストにまとめて返す
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