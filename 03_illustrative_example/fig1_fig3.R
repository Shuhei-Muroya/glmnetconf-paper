# =============================================================================
# Figure 1: solution paths of glmnet (default), glmnet (proposed) and LARS
# Figure 3: Pareto front derived from glmnet-MLP for the same dataset
#
# Data: N = 1000, p = 800, rho = 0.5, SNR = 100 (glmnetconf::data_generation)
#
# Run from the repository root.
# =============================================================================
library(glmnetconf)
library(glmnet)
library(lars)
library(MASS)
library(ggplot2)

# --- Settings (all defined here) ---
seed_value    <- 1
N_value       <- 1000
p_value       <- 800
rho_value     <- 0.5
sparse        <- 0.5
snr_value     <- 100
N_test        <- 100
t_hope        <- 20                     # time constraint for Figures 1 (b) and 3
n_coef_plot   <- 10
sample_index  <- seq_len(n_coef_plot)   # first 10 coefficients
n_s           <- 500                    # grid size for drawing the paths
s_ratio       <- 1e-4                   # smallest lambda = s_ratio * lmax
fig_width     <- 6
fig_height    <- 5
pf_fig_width  <- 9
pf_fig_height <- 6
out_dir       <- file.path("03_illustrative_example", "figures")
res_dir       <- file.path("03_illustrative_example", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# --- Data (fixed by seed) ---
set.seed(seed_value)
dat <- data_generation(N_train = N_value, N_test = N_test, p = p_value,
                       rho = rho_value, sparse_rate = sparse, snr = snr_value)
x <- dat$X_train
y <- dat$y_train
n <- nrow(x)

# --- Configuration from auto_lasso ---
result <- auto_lasso(x, y, new_x = dat$X_test, T_hope = t_hope)
saveRDS(result, file.path(res_dir, sprintf("auto_lasso_result_Thope%g.rds", t_hope)))
if (result$method != "glmnet") stop("auto_lasso selected '", result$method, "'.")
print(result$configuration[c("thresh", "nlambda")])

# =============================================================================
# Figure 1
# =============================================================================

t   <- result$configuration$thresh
lam <- result$configuration$seq_lambda

# --- Models ---
model1    <- glmnet(x = x, y = y, family = "gaussian", alpha = 1, intercept = F)  # default
lmodel    <- lars(x, y, type = "lasso")                                           # lars
model_set <- glmnet(x = x, y = y, family = "gaussian", alpha = 1,
                    lambda = lam, thresh = t, intercept = F)                      # proposed

# --- lambda grid for the coefficient paths ---
lmax  <- max(model1$lambda)
s     <- exp(seq(from = log(s_ratio * lmax), to = log(lmax), length.out = n_s))
s_log <- log(s)

# --- Coefficients ---
c_d <- t(coef(model1, s = s))
c_d <- c_d[, -1]
c_d <- c_d[, sample_index]
c_s <- t(coef(model_set, s = s))
c_s <- c_s[, -1]
c_s <- c_s[, sample_index]

#### lars
c_l <- predict.lars(lmodel, type = "coefficients", mode = "lambda", s = s * sqrt(n))$coefficients
c_l <- c_l[, sample_index]

y_min <- min(min(c_d), min(c_s), min(c_l))
y_max <- max(max(c_d), max(c_s), max(c_l))

# --- Plot and save as PDF ---
plot_path <- function(coefs, file) {
  pdf(file, width = fig_width, height = fig_height)
  matplot(s_log, coefs, type = "l",
          xlab = expression(log(lambda)), ylab = "",
          lty = 1, lwd = 3, ylim = c(y_min, y_max),
          cex.lab = 1.4, cex.axis = 1.2)
  dev.off()
}
plot_path(c_d, file.path(out_dir, "fig1a_solution_path_glmnet_default.pdf"))
plot_path(c_s, file.path(out_dir, "fig1b_solution_path_glmnet_proposed.pdf"))
plot_path(c_l, file.path(out_dir, "fig1c_solution_path_lars.pdf"))

# =============================================================================
# Figure 3
# =============================================================================
ggsave(filename = file.path(out_dir, sprintf("fig3_pareto_front_Thope%g.pdf", t_hope)),
       plot = result$Pareto_front, device = "pdf",
       width = pf_fig_width, height = pf_fig_height)
