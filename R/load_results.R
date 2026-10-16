# =============================================================================
# Loaders that flatten the chunked simulation outputs (res_chunk_*.rds) into
# data frames.
#   process_simulation_results(): Section 3.2.1 (summary dataset), Section 4.1
#   load_sensitivity_results():   Appendix D
# =============================================================================
library(dplyr)
library(purrr)

process_simulation_results <- function(directory_path) {

  files <- list.files(directory_path, pattern = "res_chunk_.*\\.rds", full.names = TRUE)
  if (length(files) == 0) {
    stop("指定されたディレクトリに該当するファイルが見つかりません。")
  }

  # --- 補助関数 ---

  # エラー記録かどうかを判定
  is_error_item <- function(item) {
    isTRUE(item$error)
  }

  # 単一の要素（1シミュレーション分）を整形する
  transform_item <- function(item) {
    lars_part <- as.data.frame(t(unlist(item$lars_list)))
    glm_part  <- item$glmnet_list$glmnet_df
    combined  <- cbind(
      lars_part[rep(1, nrow(glm_part)), , drop = FALSE],
      glm_part
    )
    list(
      lars_only = lars_part,
      combined  = combined
    )
  }

  # --- メイン処理 ---
  cat("Processing", length(files), "files...\n")

  error_log <- list()  # エラー記録を収集

  all_data_list <- map(files, function(f) {
    chunk_data <- readRDS(f)
    map(chunk_data, function(item) {
      if (is_error_item(item)) {
        # エラー記録はNULLを返してスキップ
        error_log[[length(error_log) + 1]] <<- item
        return(NULL)
      }
      tryCatch(
        transform_item(item),
        error = function(e) {
          error_log[[length(error_log) + 1]] <<- list(
            error     = TRUE,
            error_msg = conditionMessage(e),
            file      = f
          )
          NULL
        }
      )
    })
  }) %>%
    flatten() %>%
    compact()  # NULLを除去

  # エラー件数の報告
  n_errors <- length(error_log)
  if (n_errors > 0) {
    cat(sprintf("Skipped %d error item(s).\n", n_errors))
  }

  cat("Valid items:", length(all_data_list), "\n")

  # LARS単体データの行結合
  cat("Combining LARS summary...\n")
  lars_summary_matrix <- map_dfr(all_data_list, ~ .x$lars_only)
  rownames(lars_summary_matrix) <- NULL

  # LARS + glmnet 結合データの行結合
  cat("Combining All results (LARS + glmnet)...\n")
  all_combined_matrix <- map_dfr(all_data_list, ~ .x$combined)
  rownames(all_combined_matrix) <- NULL

  list(
    lars_summary = lars_summary_matrix,
    all_combined = all_combined_matrix,
    error_log    = error_log          # デバッグ用にエラー記録も返す
  )
}

load_sensitivity_results <- function(result_dir, chunk_regex) {
  files <- list.files(result_dir, pattern = chunk_regex, full.names = TRUE)
  files <- files[order(as.integer(gsub("\\D", "", basename(files))))]
  res_all  <- unlist(lapply(files, readRDS), recursive = FALSE)
  is_error <- vapply(res_all, function(x) isTRUE(x$error), logical(1))

  flatten_one <- function(x) {
    meta <- x$lars_list
    meta <- meta[vapply(meta, function(v) is.atomic(v) && length(v) == 1, logical(1))]
    df   <- x$glmnet_list$glmnet_df
    meta <- meta[setdiff(names(meta), names(df))]  # T_hope, K の重複を除く
    cbind(as.data.frame(meta, stringsAsFactors = FALSE), df)
  }

  error_df <- if (any(is_error)) {
    bind_rows(lapply(res_all[is_error], function(x) {
      data.frame(x$pattern, chunk = x$chunk, ln = x$ln,
                 error_msg = x$error_msg, stringsAsFactors = FALSE)
    }))
  } else NULL

  list(all_combined = bind_rows(lapply(res_all[!is_error], flatten_one)),
       errors = error_df,
       n_files = length(files), n_total = length(res_all), n_error = sum(is_error))
}
