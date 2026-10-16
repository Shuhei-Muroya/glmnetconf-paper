# =============================================================================
# Appendix A: effect of the configuration (tau, n_lambda) on SPE and time
#   AR(1), multivariate normal, N, p in {100, 500, 1000, 1500} with N != p,
#   100 datasets per (N, p) and case; 5 x 5 grid of configurations:
#   tau: 5 values log-equally spaced in [1e-9, 1e-7]
#   n_lambda: 5 values equally spaced in [100, 2p]
#   The paper reports case "main" (rho = 0.9, SNR = 500).
# Defines make_config_grid() and sensitivity_recipe_rep.
# =============================================================================
library(dplyr)

# ============================================================
# Common definitions
# ============================================================
n_grid_pts    <- 5
thresh_range  <- c(1e-9, 1e-7)    # log scale; the upper end is the default
nlam_min      <- 100              # default length
nlam_max_mult <- 2                # upper end = nlam_max_mult * p
maxit_val     <- 1e7
seed_upper    <- 1000000

other_fixed <- data.frame(
  sigma_struct = "AR", sprate = NA, m = NA, block_p = NA,
  beta_type = 5, beta_struct = 2, x_dis_type = 1,   # multivariate normal
  stringsAsFactors = FALSE
)

# ---- configuration grid (depends on p) ----
make_config_grid <- function(p, thresh_range, n_grid_pts, nlam_min, nlam_max_mult) {
  thresh_vec  <- 10^seq(log10(thresh_range[1]), log10(thresh_range[2]),
                        length.out = n_grid_pts)
  nlambda_vec <- round(seq(nlam_min, nlam_max_mult * p, length.out = n_grid_pts))
  grid <- expand.grid(thresh = thresh_vec, nlambda = nlambda_vec,
                      KEEP.OUT.ATTRS = FALSE)
  grid$nlambda_over_p <- grid$nlambda / p
  grid$config_id      <- seq_len(nrow(grid))
  grid
}

# ---- recipe with one dataset per row ----
# require_size: if given, keep only (n, p) pairs where n or p equals it
make_recipe <- function(N_vector, p_vector, case_settings, other_fixed,
                        n_rep, base_seed, seed_upper,
                        pattern_id_offset = 0, require_size = NULL) {
  
  Np_grid <- expand.grid(n = N_vector, p = p_vector, stringsAsFactors = FALSE)
  Np_grid <- Np_grid[Np_grid$n != Np_grid$p, ]
  if (!is.null(require_size)) {
    Np_grid <- Np_grid[Np_grid$n == require_size | Np_grid$p == require_size, ]
  }
  rownames(Np_grid) <- NULL
  
  data_recipe <- merge(Np_grid, case_settings, by = NULL)
  data_recipe <- cbind(data_recipe, other_fixed[rep(1, nrow(data_recipe)), , drop = FALSE])
  rownames(data_recipe) <- NULL
  data_recipe$pattern_id <- pattern_id_offset + seq_len(nrow(data_recipe))
  
  recipe_rep <- data_recipe[rep(seq_len(nrow(data_recipe)), each = n_rep), ]
  recipe_rep$rep_id <- rep(seq_len(n_rep), times = nrow(data_recipe))
  rownames(recipe_rep) <- NULL
  
  set.seed(base_seed)
  recipe_rep$seed <- sample(seq_len(seed_upper), nrow(recipe_rep))
  recipe_rep
}

# ============================================================
# 1. N, p in {100, 500, 1000, 1500} (N != p), main and contrast
# ============================================================
recipe_main <- make_recipe(
  N_vector      = c(100, 500, 1000, 1500),
  p_vector      = c(100, 500, 1000, 1500),
  case_settings = data.frame(case = c("main", "contrast"),
                             rho  = c(0.9, 0.3),
                             snr  = c(500, 1),
                             stringsAsFactors = FALSE),
  other_fixed   = other_fixed,
  n_rep         = 100,
  base_seed     = 2026,
  seed_upper    = seed_upper
)

sensitivity_recipe_rep <- recipe_main
cat("datasets:", nrow(sensitivity_recipe_rep), "| configs per dataset:", n_grid_pts^2, "\n")
