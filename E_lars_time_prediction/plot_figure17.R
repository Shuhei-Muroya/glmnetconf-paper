# =============================================================================
# Appendix E.2: Figure 17, predicted (lars-MLP) vs actual computation time of
# lars. The actual time is taken from the Section 4.1 results; medians over
# the 100 replications, x_dis_type = 5, N = p excluded.
#
# Run from the repository root.
# =============================================================================
source("R/load_results.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)

# 0. 設定 -------------------------------------------------------------------
series_labels    <- c(lars = "LARS", forecast = "lars-MLP")
series_levels    <- names(series_labels)
series_colors    <- setNames(c("skyblue3", "forestgreen"), series_labels)
series_linetypes <- setNames(c("dotdash", "solid"), series_labels)

LINE_SIZE    <- 1
POINT_SIZE   <- 1.5
BASE_SIZE    <- 13
X_TEXT_ANGLE <- 45
DROP_N_EQ_P  <- TRUE
LOG_Y        <- FALSE   # TRUE で縦軸を対数に

common_theme <- theme(
  legend.position = "top",
  legend.justification = "center",
  legend.box = "horizontal",
  legend.title = element_blank(),
  legend.text  = element_text(size = 14),
  legend.key.width = unit(2.5, "lines"),
  axis.title = element_text(size = 16),
  axis.text.x = element_text(angle = X_TEXT_ANGLE, hjust = 1, vjust = 1),
  strip.text = element_text(size = 12, face = "bold"),
  panel.spacing = unit(0.8, "lines"),
  plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
)

forecast_file <- file.path("E_lars_time_prediction", "results", "forecast_lars_all.csv")
save_dir      <- file.path("E_lars_time_prediction", "figures")

# 結合キー: 両方にある条件の列 + seed（all_combined 側の名前 = forecast 側の名前）
seed_map   <- c(gen_seed = "seed")
group_cols <- c("x_dis_type", "snr", "rho", "n", "p")   # 平均を取る単位（プロットと同じ）
actual_col   <- "time_lars"
forecast_col <- "T_lars_forecast"
DROP_LARS_WARNING <- FALSE   # TRUE: has_warning_lars が TRUE の行を除く

# 1. 1行ずつ結合 → (group_cols) ごとに平均 -----------------------------------
make_time_summary <- function(all_combined, forecast_df) {
  # 共通の条件列（seed 以外）をデータから取得
  pattern_cols <- setdiff(intersect(names(all_combined), names(forecast_df)),
                          c(names(seed_map), seed_map))
  cat("join keys:", c(pattern_cols, names(seed_map)), "\n")

  # 型の違い（数値 / 文字）で結合が外れないよう、キーは文字列にそろえる
  # （NA 同士は dplyr の既定でマッチする）
  key_to_chr <- function(d, cols) d %>% mutate(across(all_of(cols), as.character))

  actual <- all_combined %>%
    { if (DROP_LARS_WARNING) filter(., !as.logical(has_warning_lars)) else . } %>%
    select(all_of(c(pattern_cols, names(seed_map), actual_col))) %>%
    key_to_chr(c(pattern_cols, names(seed_map)))

  forecast <- forecast_df %>%
    rename(!!!setNames(seed_map, names(seed_map))) %>%          # seed -> gen_seed
    select(all_of(c(pattern_cols, names(seed_map), forecast_col))) %>%
    key_to_chr(c(pattern_cols, names(seed_map)))

  # 結合できなかった行の確認
  cat("actual only  :", nrow(anti_join(actual, forecast, by = c(pattern_cols, names(seed_map)))), "\n")
  cat("forecast only:", nrow(anti_join(forecast, actual, by = c(pattern_cols, names(seed_map)))), "\n")

  joined <- inner_join(actual, forecast, by = c(pattern_cols, names(seed_map)),
                       relationship = "one-to-one")              # 重複があればエラーで止める
  cat("joined rows  :", nrow(joined), "\n")

  joined %>%
    mutate(across(all_of(c("snr", "rho", "n", "p")), as.numeric)) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(time_lars     = median(.data[[actual_col]],   na.rm = TRUE),
              time_forecast = median(.data[[forecast_col]], na.rm = TRUE),
              n_rep         = dplyr::n(),
              .groups = "drop")
}


# 2. 1組み合わせ (SNR, rho) の図 --------------------------------------------
make_time_plot <- function(df_subset, title_expr = NULL) {
  p_values <- sort(unique(df_subset$p))
  p_label_levels <- paste0("p * ' ' == ' ' * ", p_values)

  plot_df <- df_subset %>%
    pivot_longer(cols = paste0("time_", series_levels), names_to = "series",
                 names_prefix = "time_", values_to = "value") %>%
    mutate(
      series  = factor(series, levels = series_levels, labels = series_labels[series_levels]),
      p_facet = factor(paste0("p * ' ' == ' ' * ", p), levels = p_label_levels)
    )

  g <- ggplot(plot_df, aes(x = n, y = value, color = series, linetype = series)) +
    geom_line(linewidth = LINE_SIZE) +
    geom_point(size = POINT_SIZE) +
    facet_wrap(~ p_facet, labeller = label_parsed, scales = "free_y", nrow = 1) +
    theme_bw(base_size = BASE_SIZE) +
    labs(x = "N", y = "Computation Time (s)", title = title_expr) +
    scale_color_manual(values = series_colors) +
    scale_linetype_manual(values = series_linetypes) +
    common_theme
  if (LOG_Y) g <- g + scale_y_log10()
  g
}

# 3. x_dis_type ごとに SNR × rho を縦に積んで保存 -----------------------------
generate_forecast_plots <- function(time_summary, save_dir, dist_values = 5,
                                    drop_n_eq_p = DROP_N_EQ_P,
                                    width_cm = 28, row_height_cm = 6, legend_height_cm = 1.2) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  if (drop_n_eq_p) time_summary <- time_summary %>% filter(n != p)

  out <- list()
  for (dist_val in dist_values) {
    df_dist <- time_summary %>% filter(x_dis_type %in% dist_val)
    combo_list <- df_dist %>% distinct(snr, rho) %>% arrange(snr, rho)
    if (nrow(combo_list) == 0) { warning(sprintf("x_dis_type = %s のデータがありません", dist_val)); next }

    plots <- map2(combo_list$snr, combo_list$rho, function(snr_val, rho_val) {
      make_time_plot(df_dist %>% filter(snr %in% snr_val, rho %in% rho_val),
                     bquote("SNR" == .(snr_val) * "," ~~ rho == .(rho_val)))
    })

    n_rows  <- length(plots)
    heights <- grid::unit(c(legend_height_cm, rep(1, n_rows)), c("cm", rep("null", n_rows)))
    combined <- wrap_plots(c(list(guide_area()), plots), ncol = 1) +
      plot_layout(guides = "collect", axis_titles = "collect", heights = heights)

    ggsave(file.path(save_dir, sprintf("Time_forecast_dist%s_ALL.pdf", dist_val)),
           combined, width = width_cm, height = legend_height_cm + row_height_cm * n_rows,
           units = "cm", limitsize = FALSE)
    out[[paste0("dist", dist_val)]] <- combined
  }
  invisible(out)
}

# --- 実行 ---
results      <- process_simulation_results(file.path("04_simulation", "results"))
all_combined <- results$all_combined %>%
  mutate(across(all_of(c("n", "p", "rho", "snr", "time_lars")), as.numeric))
forecast_df  <- read.csv(forecast_file)
time_summary <- make_time_summary(all_combined, forecast_df)
generate_forecast_plots(time_summary, save_dir = save_dir, dist_values = 5)
file.rename(file.path(save_dir, "Time_forecast_dist5_ALL.pdf"),
            file.path(save_dir, "fig17_lars_time_prediction.pdf"))
