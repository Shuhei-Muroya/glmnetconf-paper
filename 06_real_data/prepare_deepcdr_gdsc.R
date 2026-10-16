# =============================================================================
# Section 4.3.2: build the DeepCDR / GDSC dataset for one drug
# Gene expression (CCLE, 697 genes) and GDSC IC50 are read from the DeepCDR
# repository (Liu et al., 2020); cell lines with an observed IC50 are kept.
# Output: 06_real_data/data/deepcdr_gdsc1008.csv (predictors + "target")
#
# Run from the repository root.
# =============================================================================
# ---- Settings (edit here only) ----
base_url   <- "https://raw.githubusercontent.com/kimmo1019/DeepCDR/master/data/CCLE/"
expr_file  <- "genomic_expression_561celllines_697genes_demap_features.csv"
ic50_file  <- "GDSC_IC50.csv"
drug_id    <- "1008"        # GDSC drug ID (change here)

# ---- Read directly from GitHub ----
expression <- read.csv(paste0(base_url, expr_file), row.names = 1, check.names = FALSE)
ic50       <- read.csv(paste0(base_url, ic50_file),  row.names = 1, check.names = FALSE)

# ---- Build x, y for one drug ----
drug_key <- paste0("GDSC:", drug_id)
y_all    <- as.numeric(ic50[drug_key, rownames(expression)])
observed <- is.finite(y_all)

x <- as.matrix(expression[observed, ])
y <- y_all[observed]

message("N = ", nrow(x), ", p = ncol(x) = ", ncol(x))



# ---- Settings (edit here only) ----
out_dir  <- file.path("06_real_data", "data")
csv_name <- sprintf("deepcdr_gdsc%s.csv", drug_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- One CSV: predictors + target (last column) ----
target_name <- "target"
out_df <- data.frame(x, check.names = FALSE)
out_df[[target_name]] <- y

write.csv(out_df, file.path(out_dir, csv_name), row.names = FALSE)
message("saved: ", file.path(out_dir, csv_name))
