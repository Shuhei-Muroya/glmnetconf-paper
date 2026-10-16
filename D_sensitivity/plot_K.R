# =============================================================================
# Appendix D.2: Figure 14 (SPE vs K) and Figure 15 (tuning time vs K)
#
# Run from the repository root.
# =============================================================================
library(dplyr)
library(ggplot2)
source("R/load_results.R")

# ===== 設定 =====
result_dir  <- file.path("D_sensitivity", "results", "K")
save_dir    <- file.path("D_sensitivity", "figures")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
chunk_regex <- "^res_chunk_\\d+\\.rds$"
panel_w_cm  <- 6
panel_h_cm  <- 5
setting_palette <- c("red3", "skyblue3")
missing_fill    <- "grey85"

results_K <- load_sensitivity_results(result_dir, chunk_regex)
cat("Files:", results_K$n_files, " Tasks:", results_K$n_total,
    " Errors:", results_K$n_error, "\n")
if (!is.null(results_K$errors)) print(table(results_K$errors$error_msg))

all_K <- results_K$all_combined

# ===== 2. 設定ラベル (データから作成, SNR の大きい順) =====
setting_table <- all_K %>%
  distinct(rho, snr) %>%
  arrange(desc(snr), desc(rho)) %>%
  mutate(
    setting       = paste0("snr_", snr, "_rho_", rho),
    setting_label = paste0("SNR == ", snr, " * ',' ~ rho == ", rho)
  )

K_levels <- sort(unique(all_K$K))

all_K <- all_K %>%
  left_join(setting_table, by = c("rho", "snr")) %>%
  mutate(
    setting   = factor(setting, levels = setting_table$setting),
    K_f       = factor(K, levels = K_levels),
    time_total = time_tune + time_glmnet
  )

# 確認: (n, p, setting, K) ごとの件数 (想定は各 100)
print(all_K %>% count(n, p, setting, K) %>% pull(n) %>% table())

setting_colors <- setNames(setting_palette[seq_len(nrow(setting_table))],
                           setting_table$setting)
setting_labels <- setNames(parse(text = setting_table$setting_label),
                           setting_table$setting)
K_labels <- format(K_levels, big.mark = ",", scientific = FALSE, trim = TRUE)

# ===== 3. 結果が存在しない (n, p) パネル =====
missing_panels <- expand.grid(n = sort(unique(all_K$n)),
                              p = sort(unique(all_K$p))) %>%
  anti_join(distinct(all_K, n, p), by = c("n", "p"))
print(missing_panels)

# ===== 4. プロット関数 (縦軸の変数を切り替え可能) =====
make_K_boxplot <- function(df, y_var, y_label, log_y = TRUE) {
  g <- ggplot(df, aes(x = K_f, y = .data[[y_var]], fill = setting)) +
    geom_rect(data = missing_panels,
              aes(xmin = -Inf, xmax = Inf, ymin = if (log_y) 0 else -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = missing_fill) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.size = 0.6, linewidth = 0.3) +
    facet_grid(p ~ n, labeller = label_bquote(rows = p == .(p), cols = N == .(n))) +
    scale_x_discrete(labels = K_labels) +
    scale_fill_manual(values = setting_colors, labels = setting_labels) +
    theme_bw(base_size = 12) +
    labs(x = "K (number of sampled configurations)", y = y_label, fill = "Setting") +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"))
  if (log_y) g <- g + scale_y_log10()
  g
}

n_rows <- length(unique(all_K$p))
n_cols <- length(unique(all_K$n))

# ===== 5. SPE =====
p_spe <- make_K_boxplot(all_K, "SPE_glmnet", "SPE (log scale)")
print(p_spe)
ggsave(file.path(save_dir, "fig14_spe_by_K.pdf"), p_spe,
       width = panel_w_cm * n_cols, height = panel_h_cm * n_rows,
       units = "cm", limitsize = FALSE)

# ===== 6. 集計表 =====
summary_K <- all_K %>%
  group_by(setting, n, p, K) %>%
  summarise(
    n_rep          = n(),
    SPE_mean       = mean(SPE_glmnet, na.rm = TRUE),
    SPE_sd         = sd(SPE_glmnet, na.rm = TRUE),
    SPE_median     = median(SPE_glmnet, na.rm = TRUE),
    time_tune_mean = mean(time_tune, na.rm = TRUE),
    time_glmnet_mean = mean(time_glmnet, na.rm = TRUE),
    prop_within_T  = mean(time_glmnet <= T_hope, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_K, file.path("D_sensitivity", "results", "summary_K.csv"), row.names = FALSE)
print(summary_K, n = 80)



p_tune <- make_K_boxplot(all_K, "time_tune", "Tuning time (s)", log_y = FALSE)
ggsave(file.path(save_dir, "fig15_tuning_time_by_K.pdf"), p_tune,
       width = panel_w_cm * n_cols, height = panel_h_cm * n_rows,
       units = "cm", limitsize = FALSE)
