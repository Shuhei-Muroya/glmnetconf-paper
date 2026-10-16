# =============================================================================
# Appendix A: Figure 9 (SPE ratio to the default configuration) and
#             Figure 10 (computation time, with LARS as reference)
# Mean over the 100 datasets for each (N, p); rho = 0.9, SNR = 500.
#
# Run from the repository root.
# =============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 0. Settings
# ============================================================
res_dir          <- file.path("A_configuration_effect", "results")
save_dir         <- file.path("A_configuration_effect", "figures")
thresh_range     <- c(1e-9, 1e-7)    # same as in the simulation (the upper end is the default)
nlam_min         <- 100              # same as in the simulation (default length)

snr_target       <- 500              # case shown in the paper (main)
rho_target       <- 0.9

LOG_Y_RATIO      <- TRUE
LOG_Y_ABS        <- TRUE
LINE_SIZE        <- 1
POINT_SIZE       <- 1.5
BASE_SIZE        <- 13
X_TEXT_ANGLE     <- 45
fig_w_cm         <- 28
fig_h_cm         <- 24
na_panel_fill    <- "grey85"         # colour of the N = p panels

lars_label       <- "LARS"
lars_color       <- "skyblue3"
lars_linetype    <- "dotdash"

metric_defs <- data.frame(
  col   = c("SPE_glmnet", "RMSE_pred_glmnet", "RMSE_coef_glmnet", "time_glmnet"),
  label = c("SPE", "Prediction RMSE", "Coefficient RMSE", "Computation Time (s)"),
  stringsAsFactors = FALSE
)
metric_with_lars <- c("RMSE_pred_glmnet", "RMSE_coef_glmnet", "time_glmnet")

num_cols <- c("n", "p", "rho", "snr", "pattern_id", "rep_id",
              "RMSE_pred_lars", "time_lars", "RMSE_coef_lars",
              "RMSE_pred_glmnet", "time_glmnet", "SPE_glmnet", "RMSE_coef_glmnet",
              "thresh", "nlambda", "nlambda_over_p", "is_default_call", "has_warning_set")

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Read the chunks (combine lars_list and glmnet_df)
# ============================================================
files   <- list.files(res_dir, pattern = "^res_chunk_.*\\.rds$", full.names = TRUE)
res_all <- unlist(lapply(files, readRDS), recursive = FALSE)

is_err  <- vapply(res_all, function(r) isTRUE(r$error), logical(1))
cat("tasks:", length(res_all), "| failed:", sum(is_err), "\n")

all_res <- bind_rows(lapply(res_all[!is_err], function(r) {
  cbind(as.data.frame(r$lars_list, stringsAsFactors = FALSE),
        r$glmnet_list$glmnet_df)
})) %>%
  mutate(across(any_of(num_cols), as.numeric))

# ============================================================
# 2. Keep the grid rows of the target case and compute the ratio to the default corner within each dataset
# ============================================================
grid_df <- all_res %>%
  filter(is_default_call == 0,
         abs(snr - snr_target) < 1e-8,
         abs(rho - rho_target) < 1e-8)

id_cols    <- c("pattern_id", "rep_id")
base_tau   <- max(grid_df$thresh)          # tau of the default corner (taken from the data)
tau_levels <- sort(unique(grid_df$thresh))
tau_labels <- paste0("tau = ", format(tau_levels, scientific = TRUE, digits = 3))
series_levels <- c(tau_labels, lars_label)

long_df <- grid_df %>%
  select(all_of(id_cols), n, p, thresh, nlambda, nlambda_over_p,
         all_of(metric_defs$col)) %>%
  pivot_longer(all_of(metric_defs$col), names_to = "metric", values_to = "value")

base_long <- long_df %>%
  filter(thresh == base_tau, nlambda == nlam_min) %>%
  select(all_of(id_cols), metric, base = value)

ratio_long <- long_df %>%
  left_join(base_long, by = c(id_cols, "metric")) %>%
  mutate(ratio  = value / base,
         series = factor(thresh, levels = tau_levels, labels = tau_labels),
         series = factor(series, levels = series_levels),
         metric = factor(metric, levels = metric_defs$col, labels = metric_defs$label))

# colours and line types (LARS: skyblue3, dotdash)
series_colors    <- c(setNames(hcl.colors(length(tau_levels), "Dark 3"), tau_labels),
                      setNames(lars_color, lars_label))
series_linetypes <- c(setNames(rep(c("solid", "longdash", "dashed", "twodash", "dotted"),
                                   length.out = length(tau_levels)), tau_labels),
                      setNames(lars_linetype, lars_label))

# N = p panels (shaded grey)
n_values <- sort(unique(ratio_long$n))
p_values <- sort(unique(ratio_long$p))
na_cells <- expand.grid(n = n_values, p = p_values) %>% filter(n == p)

# ============================================================
# 4. Plot settings
# ============================================================
common_theme <- theme(
  legend.position      = "top",
  legend.justification = "center",
  legend.box           = "horizontal",
  legend.title         = element_blank(),
  legend.text          = element_text(size = 14),
  legend.key.width     = unit(2.5, "lines"),
  axis.title           = element_text(size = 16),
  axis.text.x          = element_text(angle = X_TEXT_ANGLE, hjust = 1, vjust = 1),
  strip.text           = element_text(size = 12, face = "bold"),
  panel.spacing        = unit(0.8, "lines"),
  plot.title           = element_text(size = 15, face = "bold", hjust = 0.5)
)

title_expr <- bquote("SNR" == .(snr_target) * "," ~~ rho == .(rho_target))

facet_labeller <- labeller(n = function(x) paste0("N = ", x),
                           p = function(x) paste0("p = ", x))

save_pdf <- function(plot, name) {
  ggsave(file.path(save_dir, name), plot,
         width = fig_w_cm, height = fig_h_cm, units = "cm", limitsize = FALSE)
}

# ============================================================
# 5. Mean over the 100 datasets
# ============================================================

# ---- mean over the 100 datasets ----
summ_mean <- ratio_long %>%
  group_by(n, p, metric, series, thresh, nlambda, nlambda_over_p) %>%
  summarise(value_mean = mean(value, na.rm = TRUE),
            ratio_mean = mean(ratio, na.rm = TRUE),
            n_rep      = dplyr::n(),
            .groups    = "drop")
cat("reps per cell (min / max):", min(summ_mean$n_rep), "/", max(summ_mean$n_rep), "\n")

# mean of the LARS reference
lars_ref_mean <- grid_df %>%
  distinct(pattern_id, rep_id, n, p, time_lars) %>%
  transmute(n, p, time_glmnet = time_lars) %>%
  pivot_longer(time_glmnet, names_to = "metric", values_to = "lars") %>%
  mutate(metric = factor(metric, levels = metric_defs$col, labels = metric_defs$label),
         series = factor(lars_label, levels = series_levels)) %>%
  group_by(n, p, metric, series) %>%
  summarise(lars_mean = mean(lars, na.rm = TRUE), .groups = "drop")

# ---- shade the N = p panels (ymin = 0 on log axes) ----
make_na_layer <- function(log_y) {
  geom_rect(data = na_cells %>% mutate(ymin_na = if (log_y) 0 else -Inf),
            inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin_na, ymax = Inf),
            fill = na_panel_fill)
}

# ---- legend labels (tau as a Greek letter) ----
tau_exp <- round(log10(tau_levels), 2)
legend_lookup <- c(lapply(tau_exp, function(e) bquote(tau == 10^.(e))),
                   list(lars_label))
names(legend_lookup) <- series_levels
legend_label_fun <- function(brks) legend_lookup[brks]

# ---- one figure ----
plot_metric_mean <- function(metric_label, type = c("ratio", "abs")) {
  type <- match.arg(type)
  d   <- summ_mean %>% filter(metric == metric_label)
  ref <- lars_ref_mean %>% filter(metric == metric_label)

  y_col <- if (type == "ratio") "ratio_mean" else "value_mean"
  y_lab <- if (type == "ratio") paste0(metric_label, " ratio to default") else metric_label
  log_y <- if (type == "ratio") LOG_Y_RATIO else LOG_Y_ABS

  g <- ggplot(d, aes(x = nlambda_over_p, y = .data[[y_col]],
                     color = series, linetype = series)) +
    make_na_layer(log_y) +
    geom_line(linewidth = LINE_SIZE) +
    geom_point(size = POINT_SIZE)

  if (type == "ratio") {
    g <- g + geom_hline(yintercept = 1, linetype = "dashed", colour = "black")
  } else if (nrow(ref) > 0) {
    g <- g + geom_hline(data = ref,
                        aes(yintercept = lars_mean, color = series, linetype = series),
                        linewidth = LINE_SIZE)
  }

  g <- g +
    facet_grid(p ~ n, labeller = facet_labeller, scales = "free_y") +
    scale_color_manual(values = series_colors,       labels = legend_label_fun) +
    scale_linetype_manual(values = series_linetypes, labels = legend_label_fun) +
    labs(x = expression(n[lambda] / p), y = y_lab) +
    theme_bw(base_size = BASE_SIZE) + common_theme

  if (log_y) g <- g + scale_y_log10()
  g
}

# ---- Figures 9 and 10 ----
save_pdf(plot_metric_mean("SPE", "ratio"), "fig9_spe_ratio.pdf")
save_pdf(plot_metric_mean("Computation Time (s)", "abs"), "fig10_computation_time.pdf")
