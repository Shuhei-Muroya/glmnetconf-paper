# =============================================================================
# Appendix E.2: predicted computation time of lars (lars-MLP) for the datasets
# of Section 4.1, via glmnetconf::forecast_lars_cptime().
# Output: E_lars_time_prediction/results/forecast_lars_all.csv
#         (chunk files forecast_chunk_<j>.rds are intermediate)
#
# Run from the repository root.
# =============================================================================
library(doParallel)
library(foreach)
library(MASS)
library(Matrix)

source("04_simulation/recipe.R")

# ---- Settings (edit here only) ----
core_n     <- 50
chunk_size <- 100
start_chunk <- 1
out_dir    <- file.path("E_lars_time_prediction", "results")
final_file <- file.path(out_dir, "forecast_lars_all.csv")
gen_files  <- "R/data_generation.R"
export_obj <- c("gen_data_flexible", "generate_sigma_wrapper",
                "coef", "get_sparse_rates", "get_top_bottom_10")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Parallel setup ----
cl <- makeCluster(core_n)
registerDoParallel(cl)

clusterExport(cl, "gen_files")
clusterEvalQ(cl, {
  library(MASS)
  library(Matrix)
  library(glmnetconf)
  for (f in gen_files) source(f)
})
clusterExport(cl, export_obj)

# ---- Chunks ----
master_recipe <- simulation_recipe_rep
n_total       <- nrow(master_recipe)
chunk_indices <- split(seq_len(n_total), ceiling(seq_len(n_total) / chunk_size))

cat("Forecasting LARS time from chunk", start_chunk, "of", length(chunk_indices), "...\n")

# ---- Main loop ----
for (j in start_chunk:length(chunk_indices)) {
  chunk_path <- file.path(out_dir, paste0("forecast_chunk_", j, ".rds"))
  if (file.exists(chunk_path)) {            # skip finished chunks (resume-friendly)
    cat("chunk", j, "already exists, skipped\n")
    next
  }
  
  Plars_j <- master_recipe[chunk_indices[[j]], ]
  LN_j    <- nrow(Plars_j)
  
  result_j <- foreach(ln = seq_len(LN_j), .packages = c("MASS", "Matrix")) %dopar% {
    tryCatch({
      pattern <- Plars_j[ln, ]
      seed    <- as.numeric(pattern$seed)
      
      # same data as the original simulation (same pattern and seed)
      data_list_ln <- gen_data_flexible(pattern, seed = seed)
      
      st <- Sys.time()
      t_forecast <- forecast_lars_cptime(data_list_ln$x)
      elapsed    <- as.numeric(difftime(Sys.time(), st, units = "secs"))
      
      data.frame(pattern,
                 T_lars_forecast = as.numeric(t_forecast),
                 forecast_elapsed = elapsed,
                 row.names = NULL)
    }, error = function(e) {
      list(error = TRUE, error_msg = conditionMessage(e),
           chunk = j, ln = ln, pattern = Plars_j[ln, ])
    })
  }
  
  saveRDS(result_j, file = chunk_path)
  
  error_tasks <- which(sapply(result_j, function(x) isTRUE(x$error)))
  if (length(error_tasks) > 0) {
    cat(sprintf("[%s] Chunk %d: %d task(s) failed (ln = %s)\n",
                Sys.time(), j, length(error_tasks), paste(error_tasks, collapse = ", ")))
  }
  cat(sprintf("[%s] Chunk %d/%d finished (%.1f%%)\n",
              Sys.time(), j, length(chunk_indices), 100 * j / length(chunk_indices)))
}

stopCluster(cl)

# ---- Combine all chunks into one CSV ----
chunk_files <- list.files(out_dir, pattern = "^forecast_chunk_\\d+\\.rds$", full.names = TRUE)
chunk_files <- chunk_files[order(as.numeric(gsub("\\D", "", basename(chunk_files))))]

all_res <- do.call(rbind, lapply(chunk_files, function(f) {
  res <- readRDS(f)
  ok  <- Filter(function(x) is.data.frame(x), res)
  do.call(rbind, ok)
}))

write.csv(all_res, final_file, row.names = FALSE)
cat("saved:", final_file, " (", nrow(all_res), "rows )\n")