# =============================================================================
# Section 3.2.1 / Appendix B.2: design of the summary dataset
#
# Builds the list of simulation settings ("recipes"). Each row is one dataset.
#   - (N, p): 15 x 15 levels, assigned by cyclic shifts (Latin-square type,
#             Fig. 11); 15 (N, p) pairs per setting
#   - 15 covariance structures x 4 design distributions x 10 beta structures
#     x 7 SNRs
#   master_recipe        : SNR in {0.05, 0.165, 0.548, 1.813, 6}  (45,000 rows)
#     master_recipe_xdis2..4 : the same split by design distribution
#     master_recipe_xdis1    : normal design, built separately (see below)
#   master_recipe_extra  : SNR in {100, 500}                      (18,000 rows)
# 63,000 datasets x 9 configurations (tau, n_lambda) = 567,000 samples.
#
# Run from the repository root:  source("01_summary_dataset/recipe.R")
# =============================================================================
library(dplyr)
set.seed(17)


# Preparation N and p -----------------------------------------------------


N_vector <- unique(c(seq(20, 200, 40), seq(200, 1000, 200), seq(1000, 3000, 400)))
p_vector <- unique(c(seq(20, 200, 40), seq(200, 1000, 200), seq(1000, 3000, 400)))


set.seed(1)
n_len <- length(N_vector)
Np_list <- do.call(rbind, lapply(0:(n_len - 1), function(shift) {
  p_indices <- ((1:n_len + shift - 1) %% n_len) + 1
  data.frame(
    n = N_vector,
    p = p_vector[p_indices],
    set_id = shift + 1  
  )
}))
print(Np_list)


# SNR ---------------------------------------------------------------------

SNR_vector <-exp(seq(log(0.05),log(6),length=5))




# Covariance structure ----------------------------------------------------

n_param_variants <- 5

rho_vector     <- c(0,0.3,0.9)       # compound symmetry / AR(1)
sprate_vector  <- c(0.01, 0.05, 0.1)       # sparse random
m_vector       <- c(1, 2, 4)                 # scale-free (Barabasi-Albert), m
block_p_vector <- c(2, 3, 5)     # block: number of blocks

# covariance structures
sigma_defs <- list(
  "1" = data.frame(sigma_struct = "compound", rho = rho_vector, sprate = NA, m = NA, block_p = NA),
  "2" = data.frame(sigma_struct = "AR", rho = rho_vector, sprate =NA, m = NA, block_p = NA),
  "3" = data.frame(sigma_struct = "random_sp", rho = NA, sprate = sprate_vector, m = NA, block_p = NA),
  "4" = data.frame(sigma_struct = "ba_hub",    rho = NA, sprate = NA, m = m_vector, block_p = NA),
  "5" = data.frame(sigma_struct = "block",     rho = NA, sprate = NA, m = NA, block_p = block_p_vector)
)

# 15 covariance settings
sigma_master <- do.call(rbind, sigma_defs)
sigma_master$sigma_id <- 1:nrow(sigma_master)


# Settings other than (N, p) ------------------------------------------------
# 15 covariance structures x 5 beta_type x 2 beta_struct x 4 distributions x 5 SNR
other_params <- expand.grid(
  sigma_idx  = 1:nrow(sigma_master),
  beta_type  = 1:5,# sparsity of beta
  beta_struct =1:2,# structure of beta: all 1 or sign * (1 + |z|)
  x_dis_type = 1:4,
  snr        = SNR_vector,
  stringsAsFactors = FALSE
)

other_params$other_id <- 1:nrow(other_params)

set.seed(123)
other_params$assigned_set_id <- rep(sample(1:15), length.out = nrow(other_params))

master_recipe <- merge(other_params, Np_list, 
                       by.x = "assigned_set_id", by.y = "set_id")


master_recipe <- merge(master_recipe, sigma_master, by.x = "sigma_idx", by.y = "sigma_id")

master_recipe <- master_recipe[, c("n", "p", "sigma_struct", "rho", "sprate", "m", "block_p", 
                                   "beta_type", "beta_struct","x_dis_type", "snr", "assigned_set_id")]
set.seed(42) # for reproducibility of the recipe
master_recipe$seed <- sample(1:1000000, nrow(master_recipe))

cat("Total number of recipt:", nrow(master_recipe), "pattern. \n")


# xdist1_normal was generated first (May 2026) from a recipe built with
# x_dis_type = 1 only, so its seeds (set.seed(42)) were drawn over 11,250 rows.
# The same steps are reproduced here.
other_params_xdis1 <- expand.grid(
  sigma_idx  = 1:nrow(sigma_master),
  beta_type  = 1:5,
  beta_struct =1:2,
  x_dis_type = 1,
  snr        = SNR_vector,
  stringsAsFactors = FALSE
)
other_params_xdis1$other_id <- 1:nrow(other_params_xdis1)

set.seed(123)
other_params_xdis1$assigned_set_id <- rep(sample(1:15), length.out = nrow(other_params_xdis1))

master_recipe_xdis1 <- merge(other_params_xdis1, Np_list,
                             by.x = "assigned_set_id", by.y = "set_id")
master_recipe_xdis1 <- merge(master_recipe_xdis1, sigma_master, by.x = "sigma_idx", by.y = "sigma_id")
master_recipe_xdis1 <- master_recipe_xdis1[, c("n", "p", "sigma_struct", "rho", "sprate", "m", "block_p",
                                               "beta_type", "beta_struct","x_dis_type", "snr", "assigned_set_id")]
set.seed(42)
master_recipe_xdis1$seed <- sample(1:1000000, nrow(master_recipe_xdis1))

master_recipe_xdis2 <- master_recipe %>% filter(x_dis_type == 2)
master_recipe_xdis3 <- master_recipe %>% filter(x_dis_type == 3)



# forth -------------------------------------------------------------------

master_recipe_xdis4 <- master_recipe %>% filter(x_dis_type == 4)



# High SNR (100, 500) -----------------------------------------------------

SNR_vector_extra<-c(100,500)
# same combinations as above
other_params_extra <- expand.grid(
  sigma_idx   = 1:nrow(sigma_master),
  beta_type   = 1:5,
  beta_struct = 1:2,
  x_dis_type  = 1:4,
  snr         = SNR_vector_extra,
  stringsAsFactors = FALSE
)
other_params_extra$other_id <- (max(other_params$other_id) + 1):(max(other_params$other_id) + nrow(other_params_extra))

set.seed(123)  # same seed as above for the assignment of assigned_set_id
other_params_extra$assigned_set_id <- rep(sample(1:15), length.out = nrow(other_params_extra))

master_recipe_extra <- merge(other_params_extra, Np_list, 
                             by.x = "assigned_set_id", by.y = "set_id")
master_recipe_extra <- merge(master_recipe_extra, sigma_master, 
                             by.x = "sigma_idx", by.y = "sigma_id")
master_recipe_extra <- master_recipe_extra[, c("n", "p", "sigma_struct", "rho", "sprate", "m", "block_p", 
                                               "beta_type", "beta_struct", "x_dis_type", "snr", "assigned_set_id")]

set.seed(4242)  # different from seed 42 used above
master_recipe_extra$seed <- sample(1:1000000, nrow(master_recipe_extra))

cat("Total number of additional recipe:", nrow(master_recipe_extra), "pattern. \n")




























