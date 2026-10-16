# =============================================================================
# Section 4.2: application to compressed sensing (Figure 8)
#
# A 32 x 32 MNIST image (label 9) is compressed with a Gaussian random
# projection to N' = 700 and reconstructed by the lasso with a two-level
# wavelet basis Psi: glmnet (default), glmnet (proposed, T_hope = 20 s), LARS.
# Output: 05_compressed_sensing/results/reconst_{def,set,lars}_N=700.csv
#         05_compressed_sensing/results/rmse_scaled.csv  (RMSE on [0, 255])
# The images in figures/ are drawn by plot_reconstruction.py.
#
# Run from the repository root.
# =============================================================================
library(MASS)
library(glmnetconf)
library(lars)
library(glmnet)
compses_glmnetconf <- function(original,Psi,N=256,seed=1,save=FALSE,folder_name = "result"){
  set.seed(seed)
  p<-length(original)
  mu<-rep(0,p)
  sigma<-diag(1,p,p)
  Z<-mvrnorm(N,mu,sigma)
  y<-Z%*%original
  X<-Z%*%Psi
  
  # lars
  a<-Sys.time()
  cvlars_p<-lars::cv.lars(X,y,plot.it = FALSE,type = "lasso" ,max.step=p,use.Gram=FALSE,intercept=FALSE)
  lars.model_p<-lars::lars(X,y,use.Gram=FALSE,type = "lasso",trace = FALSE,intercept=FALSE)
  b<-Sys.time()
  time_lars<-b-a
  
  
  coef_p<-coef.lars(lars.model_p,mode="lambda",s=cvlars_p$index[which.min(cvlars_p$cv)])#,mode="lambda"
  theta_reconst2<-Psi%*%coef_p
  
  
  #glmnet, デフォルト
  set.seed(seed)
  a<-Sys.time()
  cv_def<-cv.glmnet(X,y,family = "gaussian",alpha = 1,intercept=FALSE)
  model_def<-glmnet(X,y,family = "gaussian",alpha = 1,intercept=FALSE)
  b<-Sys.time()
  time_def<-b-a
  coef_def<-stats::coef(model_def,s=cv_def$lambda.min)
  coef_def<-as.numeric(coef_def)[-1]#interseptは0になるFALSEで、だからここで抜く。FALSEにしても消えるわけではない
  theta_reconst_def<-Psi%*%coef_def
  
  ##auto
  lam<-numeric()
  set<-tune_configuration(X,T_hope=20,size = 1000)
  pl<-set$best_configuration$nlambda
  t<-set$best_configuration$thresh
  minn<-min(model_def$lambda)####修正
  len_model1_lambda<-length(model_def$lambda)
  model1_lambda<-model_def$lambda
  #print(model1_lambda)
  if(pl==0){
    lam<-model1_lambda
  }else if(pl<len_model1_lambda){
    add<-seq(minn,0,by=-(minn/(100-len_model1_lambda)))[-1]#100まで増やす必要ない。
    #単に長さplにする
    add<-seq(minn,0,length.out=pl)#[-1]
    lam<-add##ここの修正
  }else{
    add<-seq(minn,0,by=-(minn/(pl-len_model1_lambda)))[-1]#続きを増やしている。
    lam<-c(model1_lambda,add)#増やしたlambdaの候補
  }
  
  pl_l<-length(lam)
  st<-Sys.time()
  cv_set<-cv.glmnet(X,y,family = "gaussian",alpha = 1,intercept=FALSE,thresh=t,lambda=lam)
  model1_set<-glmnet(X,y,family = "gaussian",alpha = 1,intercept=FALSE,thresh=t,lambda=lam)
  gt<-Sys.time()
  time_set<-gt-st
  coef_set<-stats::coef(model1_set,s=cv_set$lambda.min)
  coef_set<-as.numeric(coef_set)[-1]
  theta_reconst_set<-Psi%*%coef_set
  
  if (save) {
    dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
    larsname<-paste0("reconst_lars_N=",N,".csv")
    defname<-paste0("reconst_def_N=",N,".csv")
    setname<-paste0("reconst_set_N=",N,".csv")
    write.csv(theta_reconst2, file.path(folder_name, larsname), row.names = FALSE)
    write.csv(theta_reconst_def, file.path(folder_name, defname), row.names = FALSE)
    write.csv(theta_reconst_set, file.path(folder_name, setname), row.names = FALSE)
  }
  
  result<-list(theta_lars=theta_reconst2,
               theta_def=theta_reconst_def,
               theta_set=theta_reconst_set,
               time_lars=time_lars,
               time_def=time_def,
               time_set=time_set)
  return(result)
}

# 実行 ----------------------------------------------------------------------
data_dir    <- file.path("05_compressed_sensing", "data")
results_dir <- file.path("05_compressed_sensing", "results")

Psi<-read.csv(file.path(data_dir, "waveletbasis.csv"))
Psi<-as.matrix(Psi)
Psi<-t(Psi)

file_path <- file.path(data_dir, "mnist_label_9_image_4.csv")
data <- read.csv(file_path, header = FALSE)
original<-as.numeric(unlist(t(data)))#tはplotとかと合わせるために必要

a<-compses_glmnetconf(original,Psi,N=700,save=T,folder_name = results_dir)
print(a[c("time_lars", "time_def", "time_set")])


# RMSEまで計算 ----------------------------------------------------------------

library(dplyr)
original<-read.csv(file.path(data_dir, "mnist_label_9_image_4.csv"),header = F)
def<-read.csv(file.path(results_dir, "reconst_def_N=700.csv"))
lars<-read.csv(file.path(results_dir, "reconst_lars_N=700.csv"))
set<-read.csv(file.path(results_dir, "reconst_set_N=700.csv"))




# スケーリング ------------------------------------------------------------------


# 1. まずスケーリング関数を定義
rescale_255 <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (x_max == x_min) return(numeric(length(x)))
  
  # 0-255に変換
  return((x - x_min) / (x_max - x_min) * 255)
}

# 2. original と各データをベクトル化
original_vec <- as.vector(t(original))
results_raw <- list(
  def  = unlist(def),
  lars = unlist(lars),
  set  = unlist(set)
)

# 3. すべてのデータを個別にスケーリング
# originalも比較対象としてスケーリングしておく必要があります
original_scaled <- rescale_255(original_vec)
results_scaled <- lapply(results_raw, rescale_255)

# 4. RMSEを計算する関数の定義
calc_rmse <- function(target, ref) {
  sqrt(mean((target - ref)^2, na.rm = TRUE))
}

# 5. スケーリング済みデータ同士でRMSEを一括計算
rmses_scaled <- sapply(results_scaled, calc_rmse, ref = original_scaled)

# 結果の表示
print("--- Scaled RMSE (0-255) ---")
print(rmses_scaled)


# 結果の保存 -------------------------------------------------------------------

# 1. 保存先のパスを指定
save_path <- file.path(results_dir, "rmse_scaled.csv")

# 2. 結果をデータフレームに整形（手法名とRMSE値）
rmse_df <- data.frame(
  Method = names(rmses_scaled),
  RMSE_Scaled = as.numeric(rmses_scaled)
)

# 3. 保存
write.csv(rmse_df, save_path, row.names = FALSE)

# 確認表示
cat("--- 以下のパスにRMSEの結果を保存しました ---\n")
print(save_path)
print(rmse_df)


