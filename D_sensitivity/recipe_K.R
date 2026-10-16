# =============================================================================
# Appendix D.2: recipe for the sensitivity to K (Figures 14, 15)
#   Same as recipe_T_hope.R but T_hope = 10 fixed and K in {100, 1000, 10000}.
#   Defines sensitivity_recipe_rep.
# =============================================================================
library(dplyr)


N_vector <- c(100, 500, 1000, 1500)
p_vector <- c(100, 500, 1000, 1500)

Np_grid <- expand.grid(n = N_vector, p = p_vector, stringsAsFactors = FALSE)
Np_grid <- Np_grid[Np_grid$n != Np_grid$p, ]  # exclude N = p
rownames(Np_grid) <- NULL

cat("Number of (n, p) combinations (N != p):", nrow(Np_grid), "\n")
print(Np_grid)

T_hope_fixed <- 10
K_vector     <- c(100, 1000, 10000)

case_settings <- data.frame(
  case = c("main", "contrast"),
  rho  = c(0.9, 0.3),
  snr  = c(500, 1),
  stringsAsFactors = FALSE
)

other_fixed <- data.frame(
  sigma_struct = "AR",
  sprate       = NA,
  m            = NA,
  block_p      = NA,
  beta_type    = 5,
  beta_struct  = 2,
  x_dis_type   = 5,
  stringsAsFactors = FALSE
)

# grid of K x case (T_hope fixed)
k_case_grid <- expand.grid(
  K      = K_vector,
  case   = case_settings$case,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(T_hope = T_hope_fixed) %>%
  dplyr::left_join(case_settings, by = "case")

# cross join Np_grid x (K x case)
sensitivity_recipe <- merge(Np_grid, k_case_grid, by = NULL)
sensitivity_recipe <- cbind(sensitivity_recipe, other_fixed[rep(1, nrow(sensitivity_recipe)), , drop = FALSE])
rownames(sensitivity_recipe) <- NULL

sensitivity_recipe$pattern_id <- seq_len(nrow(sensitivity_recipe))

n_rep <- 100

sensitivity_recipe_rep <- sensitivity_recipe[rep(sensitivity_recipe$pattern_id, each = n_rep), ]
sensitivity_recipe_rep$rep_id <- rep(seq_len(n_rep), times = nrow(sensitivity_recipe))
rownames(sensitivity_recipe_rep) <- NULL

set.seed(2026)
sensitivity_recipe_rep$seed <- sample(1:1000000, nrow(sensitivity_recipe_rep))

cat("Number of (n, p, case, K) combinations:", nrow(sensitivity_recipe), "\n")
cat("Reps per combination:", n_rep, "\n")
cat("Total number of simulations:", nrow(sensitivity_recipe_rep), "\n")

print(sensitivity_recipe %>% dplyr::select(n, p, case, rho, snr, T_hope, K) %>% head(20))
head(sensitivity_recipe_rep)
