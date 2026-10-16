# =============================================================================
# Section 4.1: simulation settings
#   N, p in {100, 1000, 1500, 2000, 2500}, AR(1) with rho in {0.3, 0.9},
#   SNR in {1, 100, 500}, floor(p/2) nonzero coefficients sign * (1 + |z|),
#   100 replications per setting.
#   x_dis_type = 1:5 are all simulated; the paper reports x_dis_type = 5
#   (multivariate t with 5 degrees of freedom, not used in the summary dataset).
# Defines simulation_recipe_rep. Also used by E_lars_time_prediction/.
# =============================================================================
# Preparation N and p -----------------------------------------------------
N_vector <- c(100, 1000, 1500,2000,2500)
p_vector <- c(100, 1000, 1500,2000,2500)
Np_list <- expand.grid(n = N_vector, p = p_vector, stringsAsFactors = FALSE)
print(Np_list)

# SNR ---------------------------------------------------------------------
#SNR_vector <- c(6)

SNR_vector <- c(1,100, 500)#  high snr

# Covariance structure ----------------------------------------------------
rho_vector <- c(0.3, 0.9)       # AR(1) 用

sigma_defs <- list(
  "1" = data.frame(sigma_struct = "AR", rho = rho_vector, sprate = NA, m = NA, block_p = NA)
)
sigma_master <- do.call(rbind, sigma_defs)
sigma_master$sigma_id <- 1:nrow(sigma_master)

# Total simulation 数 ------------------------------------------------------
other_params <- expand.grid(
  sigma_idx   = 1:nrow(sigma_master),
  beta_type   = 5,   # betaのスパース割合, 半分がsparseの状況にする。
  beta_struct = 2,     # betaの構造, 符号付+乱数込みにする
  x_dis_type  = 1:5,
  snr         = SNR_vector,
  stringsAsFactors = FALSE
)
other_params$other_id <- 1:nrow(other_params)

# Np_list × other_params の完全クロス結合(全組み合わせを網羅)
master_recipe <- merge(other_params, Np_list, by = NULL)  # by=NULLでcross join
master_recipe <- merge(master_recipe, sigma_master, by.x = "sigma_idx", by.y = "sigma_id")

master_recipe <- master_recipe[, c("n", "p", "sigma_struct", "rho", "sprate", "m", "block_p",
                                   "beta_type", "beta_struct", "x_dis_type", "snr")]

set.seed(42) # レシピ生成自体の再現性のため
master_recipe$seed <- sample(1:1000000, nrow(master_recipe))

cat("Total number of recipe:", nrow(master_recipe), "pattern. \n")

simulation_recipe <- master_recipe



# 各パターンをn_rep回繰り返すための設定 -------------------------------------
n_rep <- 100  # 繰り返し回数

# 元のパターンにIDを付与(後で集計時にグループ化するため)
simulation_recipe$pattern_id <- seq_len(nrow(simulation_recipe))

# 各行をn_rep回複製
simulation_recipe_rep <- simulation_recipe[rep(simulation_recipe$pattern_id, each = n_rep), ]
simulation_recipe_rep$rep_id <- rep(seq_len(n_rep), times = nrow(simulation_recipe))
rownames(simulation_recipe_rep) <- NULL

# 複製ごとに新しいseedを付与(元のseed列は上書き)
set.seed(2026)  # seed生成自体の再現性のため
simulation_recipe_rep$seed <- sample(1:1000000, nrow(simulation_recipe_rep))

# 確認
cat("Number of patterns:", nrow(simulation_recipe), "\n")
cat("Reps per pattern:", n_rep, "\n")
cat("Total number of simulations:", nrow(simulation_recipe_rep), "\n")

head(simulation_recipe_rep)
