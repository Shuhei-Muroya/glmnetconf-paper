# =============================================================================
# Section 4.1: Figures 4-7
#   Fig. 4 prediction RMSE, Fig. 5 SPE, Fig. 6 computation time,
#   Fig. 7 RMSE of the coefficients
# Only x_dis_type = 5 (multivariate t) is plotted; N = p is excluded.
#
# Run from the repository root.
# =============================================================================
source("R/load_results.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)

# 0. Settings -----------------------------------------------------------------
method_labels <- c(
  default = "glmnet (default)",
  glmnet  = "glmnet (proposed)",
  lars    = "LARS"
)
method_levels <- names(method_labels)
method_colors <- c("gray40", "red3", "skyblue3")
names(method_colors) <- method_labels[method_levels]

# per metric: column prefix, y-axis label and output file name;
# methods are the columns "prefix + method_levels" present in the data
# (SPE has no LARS column, so only default and proposed are drawn)
metric_specs <- list(
  RMSE = list(prefix = "RMSE_pred_", ylab = "RMSE",                 file = "RMSE_comparison"),
  Time = list(prefix = "time_",      ylab = "Computation Time (s)", file = "Time_comparison"),
  Coef = list(prefix = "RMSE_coef_", ylab = "RMSE of Coefficients", file = "Coef_RMSE_comparison"),
  SPE  = list(prefix = "SPE_",       ylab = "SPE",                  file = "SPE_comparison")
)

LINE_SIZE  <- 1
POINT_SIZE <- 1.5
BASE_SIZE  <- 13
X_TEXT_ANGLE <- 45   # angle of the x-axis (N) tick labels; 0 = horizontal

# points excluded from the plots
DROP_N_EQ_P <- TRUE  # TRUE: drop N = p / FALSE: plot all

common_theme <- theme(
  legend.position = "top",
  legend.justification = "center",
  legend.box = "horizontal",
  legend.title = element_text(size = 14, face = "bold"),
  legend.text  = element_text(size = 14),
  legend.key.width = unit(2.5, "lines"),
  axis.title = element_text(size = 16),
  axis.text.x = element_text(angle = X_TEXT_ANGLE, hjust = 1, vjust = 1),
  strip.text = element_text(size = 12, face = "bold"),
  panel.spacing = unit(0.8, "lines"),
  plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
)

# 1. Plot of one metric for one (SNR, rho) combination ---------------------------
make_metric_plot <- function(df_subset, spec, title_expr = NULL) {
  value_cols <- intersect(paste0(spec$prefix, method_levels), names(df_subset))

  p_values <- sort(unique(df_subset$p))
  p_label_levels <- paste0("p * ' ' == ' ' * ", p_values)

  plot_df <- df_subset %>%
    dplyr::group_by(n, p) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(value_cols), ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(value_cols), names_to = "method",
                        names_prefix = spec$prefix, values_to = "value") %>%
    dplyr::mutate(
      method  = factor(method, levels = method_levels, labels = method_labels[method_levels]),
      p_facet = factor(paste0("p * ' ' == ' ' * ", p), levels = p_label_levels)
    )

  ggplot(plot_df, aes(x = n, y = value, color = method, linetype = method)) +
    geom_line(linewidth = LINE_SIZE) +
    geom_point(size = POINT_SIZE) +
    facet_wrap(~ p_facet, labeller = label_parsed, scales = "free_y", nrow = 1) +
    theme_bw(base_size = BASE_SIZE) +
    labs(x = "N", y = spec$ylab, color = "Method", linetype = "Method", title = title_expr) +
    scale_color_manual(values = method_colors) +
    common_theme
}

# 2. For each x_dis_type, stack the (SNR, rho) rows into one figure (_ALL) and save it ----
generate_all_plots <- function(all_combined, save_dir, dist_values = 5,
                               drop_n_eq_p = DROP_N_EQ_P,
                               width_cm = 28, row_height_cm = 6, legend_height_cm = 1.2) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  if (drop_n_eq_p) {
    n_before <- nrow(all_combined)
    all_combined <- all_combined %>% dplyr::filter(n != p)
    cat(sprintf("Dropped %d rows with N = p\n", n_before - nrow(all_combined)))
  }

  combined_plots <- list()

  for (dist_val in dist_values) {
    df_dist <- all_combined %>% dplyr::filter(x_dis_type %in% dist_val)

    # (SNR, rho) combinations present in the data
    combo_list <- df_dist %>%
      dplyr::distinct(snr, rho) %>%
      dplyr::arrange(snr, rho)
    print(combo_list)

    if (nrow(combo_list) == 0) {
      warning(sprintf("no data for x_dis_type = %s", dist_val))
      next
    }

    # subsets and titles for each combination (rho as a Greek letter)
    subsets <- purrr::map2(combo_list$snr, combo_list$rho, function(snr_val, rho_val) {
      df_dist %>% dplyr::filter(snr %in% snr_val, rho %in% rho_val)   # %in% also handles NA
    })
    titles <- purrr::map2(combo_list$snr, combo_list$rho, function(snr_val, rho_val) {
      bquote("SNR" == .(snr_val) * "," ~~ rho == .(rho_val))
    })

    for (metric in names(metric_specs)) {
      spec <- metric_specs[[metric]]

      plots <- purrr::map2(subsets, titles, function(d, ttl) make_metric_plot(d, spec, ttl))

      # legend row (guide_area) on top, then one row per (SNR, rho):
      #   legend -> SNR title -> panels
      # guides = "collect"      : one shared legend in guide_area
      # axis_titles = "collect" : one y-axis title on the left and one x-axis title (N) at the bottom
      #                           (requires patchwork >= 1.2.0)
      n_rows  <- length(plots)
      heights <- grid::unit(c(legend_height_cm, rep(1, n_rows)),
                            c("cm", rep("null", n_rows)))

      combined <- patchwork::wrap_plots(c(list(patchwork::guide_area()), plots), ncol = 1) +
        patchwork::plot_layout(guides = "collect", axis_titles = "collect", heights = heights)

      ggsave(file.path(save_dir, sprintf("%s_dist%s_ALL.pdf", spec$file, dist_val)),
             combined, width = width_cm, height = legend_height_cm + row_height_cm * n_rows,
             units = "cm", limitsize = FALSE)

      combined_plots[[sprintf("dist%s_%s", dist_val, metric)]] <- combined
    }
  }

  invisible(combined_plots)
}

# 3. Run -----------------------------------------------------------------------
results      <- process_simulation_results(file.path("04_simulation", "results"))
all_combined <- results$all_combined

cols_to_numeric <- c(
  "n", "p", "rho", "sprate", "m", "block_p", "snr", "error_sd",
  "RMSE_pred_default", "time_default",
  "RMSE_pred_lars", "time_lars",
  "SPE_default", "sd_SPE_default",
  "RMSE_coef_default", "RMSE_coef_lars",
  "nlambda", "thresh",
  "RMSE_pred_glmnet", "time_glmnet",
  "SPE_glmnet", "sd_SPE_glmnet",
  "RMSE_coef_glmnet", "best_lambda"
)
all_combined <- all_combined %>%
  mutate(across(all_of(cols_to_numeric), as.numeric))

save_dir <- file.path("04_simulation", "figures")
results_plots <- generate_all_plots(all_combined, save_dir = save_dir, dist_values = 5)

# rename to the figure numbers of the paper
fig_names <- c(RMSE_comparison      = "fig4_prediction_rmse",
               SPE_comparison       = "fig5_spe",
               Time_comparison      = "fig6_computation_time",
               Coef_RMSE_comparison = "fig7_coefficient_rmse")
for (k in names(fig_names)) {
  file.rename(file.path(save_dir, sprintf("%s_dist5_ALL.pdf", k)),
              file.path(save_dir, sprintf("%s.pdf", fig_names[[k]])))
}
