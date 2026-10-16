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

SNR_vector <- c(1,100, 500)

# Covariance structure ----------------------------------------------------
rho_vector <- c(0.3, 0.9)       # AR(1)

sigma_defs <- list(
  "1" = data.frame(sigma_struct = "AR", rho = rho_vector, sprate = NA, m = NA, block_p = NA)
)
sigma_master <- do.call(rbind, sigma_defs)
sigma_master$sigma_id <- 1:nrow(sigma_master)

# Settings other than (N, p) ------------------------------------------------
other_params <- expand.grid(
  sigma_idx   = 1:nrow(sigma_master),
  beta_type   = 5,   # floor(0.5p) nonzero coefficients
  beta_struct = 2,     # nonzero coefficients sign * (1 + |z|)
  x_dis_type  = 1:5,
  snr         = SNR_vector,
  stringsAsFactors = FALSE
)
other_params$other_id <- 1:nrow(other_params)

# all combinations of Np_list x other_params
master_recipe <- merge(other_params, Np_list, by = NULL)  # cross join
master_recipe <- merge(master_recipe, sigma_master, by.x = "sigma_idx", by.y = "sigma_id")

master_recipe <- master_recipe[, c("n", "p", "sigma_struct", "rho", "sprate", "m", "block_p",
                                   "beta_type", "beta_struct", "x_dis_type", "snr")]

set.seed(42) # for reproducibility of the recipe
master_recipe$seed <- sample(1:1000000, nrow(master_recipe))

cat("Total number of recipe:", nrow(master_recipe), "pattern. \n")

simulation_recipe <- master_recipe



# Replications --------------------------------------------------------------
n_rep <- 100  # number of replications

# pattern ID (used to group replications)
simulation_recipe$pattern_id <- seq_len(nrow(simulation_recipe))

# replicate each row n_rep times
simulation_recipe_rep <- simulation_recipe[rep(simulation_recipe$pattern_id, each = n_rep), ]
simulation_recipe_rep$rep_id <- rep(seq_len(n_rep), times = nrow(simulation_recipe))
rownames(simulation_recipe_rep) <- NULL

# new seed for each replication (overwrites the seed column)
set.seed(2026)  # for reproducibility of the seeds
simulation_recipe_rep$seed <- sample(1:1000000, nrow(simulation_recipe_rep))

cat("Number of patterns:", nrow(simulation_recipe), "\n")
cat("Reps per pattern:", n_rep, "\n")
cat("Total number of simulations:", nrow(simulation_recipe_rep), "\n")

head(simulation_recipe_rep)
