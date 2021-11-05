#' @export
get_pipeline_df <- function(dir_base, dir_sub = "", fn, df_cols = NULL, obj_function) {
  dir_base_sub <- paste0(dir_base, "/", paste0(dir_sub, collapse = "/"))
  file_vec <- list.files(dir_base_sub)
  obj_ind <- fn %in% file_vec

  if (!obj_ind) {
    dir_vec <- list.dirs(dir_base_sub,
      full.names = FALSE,
      recursive = FALSE
    )
    if (any(dir_vec == "mod_full")) dir_vec <- "mod_full"
    if (any(dir_vec == "results")) dir_vec <- "results"
    if (length(dir_vec) == "0") {
      if (is.null(df_cols)) {
        return(NULL)
      }
      out_mat <- matrix(rep(NA, length(df_cols)), ncol = length(df_cols))
      colnames(out_mat) <- df_cols
      out_df <- out_mat %>% as_tibble()
      return(out_df)
    }
    if (identical(dir_sub, "")) {
      dir_sub_list <- map(dir_vec, function(x) x)
    } else {
      dir_sub_list <- map(dir_vec, function(x) c(dir_sub, x))
    }
    out_list <- map_df(dir_sub_list, function(x) {
      get_pipeline_df(
        dir_base = dir_base, dir_sub = x,
        fn = fn, df_cols = df_cols
      )
    })
  } else {
    obj <- readRDS(file.path(dir_base_sub, fn))
    mat_var <- matrix(rep(dir_sub, nrow(obj)),
      byrow = TRUE,
      ncol = length(dir_sub)
    )
    colnames(mat_var) <- paste0("V", seq.int(1, length(dir_sub)))
    df_var <- as_tibble(mat_var)
    out_list <- bind_cols(df_var, obj)
  }

  out_list
}
