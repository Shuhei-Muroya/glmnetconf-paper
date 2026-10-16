# =============================================================================
# Appendix D.1: Figure 12 (glmnet computation time vs T_hope) and
#               Figure 13 (SPE vs T_hope)
#
# Run from the repository root.
# =============================================================================
library(dplyr)
library(ggplot2)
source("R/load_results.R")

# ===== settings =====
result_dir  <- file.path("D_sensitivity", "results", "T_hope")
save_dir    <- file.path("D_sensitivity", "figures")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

chunk_regex <- "^res_chunk_\\d+\\.rds$"
panel_w_cm  <- 6
panel_h_cm  <- 5
setting_palette <- c("red3", "skyblue3")
missing_fill    <- "grey85"
ref_shape       <- 4     # marker of the specified T_hope (x)
ref_size        <- 2
ref_stroke      <- 0.8

results_Thope <- load_sensitivity_results(result_dir, chunk_regex)
cat("Files:", results_Thope$n_files, " Tasks:", results_Thope$n_total,
    " Errors:", results_Thope$n_error, "\n")
if (!is.null(results_Thope$errors)) print(table(results_Thope$errors$error_msg))

all_Thope <- results_Thope$all_combined

# ===== 2. Setting labels (from the data, larger SNR first) =====
setting_table <- all_Thope %>%
  distinct(rho, snr) %>%
  arrange(desc(snr), desc(rho)) %>%
  mutate(
    setting       = paste0("snr_", snr, "_rho_", rho),
    setting_label = paste0("SNR == ", snr, " * ',' ~ rho == ", rho)
  )

T_hope_levels <- sort(unique(all_Thope$T_hope))

all_Thope <- all_Thope %>%
  left_join(setting_table, by = c("rho", "snr")) %>%
  mutate(
    setting  = factor(setting, levels = setting_table$setting),
    T_hope_f = factor(T_hope, levels = T_hope_levels)
  )

# check: number of rows per (n, p, setting, T_hope) (100 each)
print(all_Thope %>% count(n, p, setting, T_hope) %>% pull(n) %>% table())

setting_colors <- setNames(setting_palette[seq_len(nrow(setting_table))],
                           setting_table$setting)
setting_labels <- setNames(parse(text = setting_table$setting_label),
                           setting_table$setting)
T_hope_labels  <- format(T_hope_levels, big.mark = ",", scientific = FALSE, trim = TRUE)

# ===== 3. (n, p) panels without results, and reference points =====
missing_panels <- expand.grid(n = sort(unique(all_Thope$n)),
                              p = sort(unique(all_Thope$p))) %>%
  anti_join(distinct(all_Thope, n, p), by = c("n", "p"))
print(missing_panels)

ref_df <- all_Thope %>% distinct(n, p, T_hope, T_hope_f)  # panels with results only

# ===== 4. Plot function (x and y variables selectable) =====
make_sensi_boxplot <- function(df, x_var, x_labels, x_title,
                               y_var, y_label, log_y = TRUE, ref = NULL) {
  g <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = setting)) +
    geom_rect(data = missing_panels,
              aes(xmin = -Inf, xmax = Inf, ymin = if (log_y) 0 else -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = missing_fill) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.size = 0.6, linewidth = 0.3)
  if (!is.null(ref)) {
    g <- g + geom_point(data = ref, aes(x = .data[[x_var]], y = T_hope),
                        inherit.aes = FALSE,
                        shape = ref_shape, size = ref_size, stroke = ref_stroke)
  }
  g <- g +
    facet_grid(p ~ n, labeller = label_bquote(rows = p == .(p), cols = N == .(n))) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_manual(values = setting_colors, labels = setting_labels) +
    theme_bw(base_size = 12) +
    labs(x = x_title, y = y_label, fill = "Setting") +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"))
  if (log_y) g <- g + scale_y_log10()
  g
}

n_rows <- length(unique(all_Thope$p))
n_cols <- length(unique(all_Thope$n))
x_title_Thope <- expression(T[hope] ~ "(s)")

# ===== 5. Computation time of glmnet (x = T_hope) =====
p_time <- make_sensi_boxplot(all_Thope, "T_hope_f", T_hope_labels, x_title_Thope,
                             "time_glmnet", "Computation time (s, log scale)",
                             log_y = TRUE, ref = ref_df)
print(p_time)
ggsave(file.path(save_dir, "fig12_time_by_Thope.pdf"), p_time,
       width = panel_w_cm * n_cols, height = panel_h_cm * n_rows,
       units = "cm", limitsize = FALSE)

# ===== 6. SPE =====
p_spe <- make_sensi_boxplot(all_Thope, "T_hope_f", T_hope_labels, x_title_Thope,
                            "SPE_glmnet", "SPE (log scale)", log_y = TRUE)
print(p_spe)
ggsave(file.path(save_dir, "fig13_spe_by_Thope.pdf"), p_spe,
       width = panel_w_cm * n_cols, height = panel_h_cm * n_rows,
       units = "cm", limitsize = FALSE)

# ===== 7. Summary table =====
summary_Thope <- all_Thope %>%
  group_by(setting, n, p, T_hope) %>%
  summarise(
    n_rep            = n(),
    SPE_mean         = mean(SPE_glmnet, na.rm = TRUE),
    SPE_sd           = sd(SPE_glmnet, na.rm = TRUE),
    SPE_median       = median(SPE_glmnet, na.rm = TRUE),
    time_glmnet_mean = mean(time_glmnet, na.rm = TRUE),
    time_ratio_mean  = mean(time_glmnet / T_hope, na.rm = TRUE),
    prop_within_T    = mean(time_glmnet <= T_hope, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_Thope, n = 100)
