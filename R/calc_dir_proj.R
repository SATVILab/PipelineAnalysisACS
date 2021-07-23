#' @title Get project directory based on iterators
#'
#' @export
calc_dir_proj <- function(dir_base, iter = list(), ...) {
  d_list <- iter %>% append(list(...))
  dir_sub <- purrr::map_chr(
    seq_along(d_list),
    function(i) {
      get_folder(nm = names(d_list)[i],
                 value = d_list[[i]])
    }
  ) %>%
    paste0(collapse = "/")
  file.path(dir_base, dir_sub) %>%
    fs::path_norm()
}

get_folder <- function(nm, value) {
  switch(nm,
    "max_ttb" = as.character(value),
    "conf" = ifelse(
      identical(value, "none"), "no_conf", paste0(value, collapse = "")
      ),
    "rem_il6" = ifelse(value, "il6_exc", "il6_inc"),
    "wins" = ifelse(value, "wins", "wins_n"),
    "equi_d_knots" = ifelse(value, "equi_d", "equi_d_n"),
    "ds" =,
    "dataset" = value,
    "stim" = value,
    "cyt_response_type" = value,
    "cyt_response_type_grp" = value,
    stop(paste0("name ", nm, "not recognised in calc_dir_proj"))
  )
}
