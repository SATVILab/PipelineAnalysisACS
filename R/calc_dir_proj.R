#' @title Get project directory based on iterators
#'
#' @export
calc_dir_proj <- function(dir_base, iter = list(), ...) {
  d_list <- iter %>% append(list(...))
  dir_sub <- purrr::map_chr(
    seq_along(d_list),
    function(i) {
      out <- get_folder(nm = names(d_list)[i],
                 value = d_list[[i]])
      out
    }
  ) %>%
    paste0(collapse = "/")
  file.path(dir_base, dir_sub) %>%
    fs::path_norm()
}

get_folder <- function(nm, value) {
  value <- switch(as.character(is.na(value)),
                  "TRUE" = "-",
                  "FALSE" = value)
  out <- switch(nm,
    "ttb_max" = as.character(value),
    "var_conf" = ,
    "conf" = ifelse(
      identical(value[[1]], "none"), "no_conf", paste0(
        stringr::str_sub(value[[1]], end = min(nchar(value[[1]]), 3)),
                         collapse = "")
      ),
    "pkg" = value,
    "rem_il6" = ifelse(value, "il6_e", "il6_i"),
    "wins" = ifelse(value, "wins", "wins_n"),
    "equi_d_knots" = ifelse(value, "eq_d", "eq_d_n"),
    "ds" =,
    "dataset" = value,
    "stim" = value,
    "var_dep" = value,
    "cyt_response_type" = value,
    "var_offset" = paste0("o-", stringr::str_sub(value, end = min(nchar(value), 3)) %>%
                            paste0(collapse = "_")),
    "cyt_response_type_grp" = value,
    "var_exp" = paste0("exp-", paste0(
      stringr::str_sub(value, end = min(nchar(value), 5)),
      collapse = "")),
    "var_int" = ifelse(
      is.null(value[[1]]),
      "int_n",
      "int"),
    "var_exp_spline" =
      switch(
        as.character(identical(value, "none")),
        "TRUE" = "exp_s-none",
        "FALSE" = paste0(

          "exp_s-",
          paste0(purrr::map_chr(seq_along(value), function(i) {
            nm_curr <- names(value)[i]
            if(grepl("^tc~", nm_curr)) {
              nm_curr <- stringr::str_split(nm_curr, "~")[[1]][3]
            }
            out <- stringr::str_sub(nm_curr, end = min(6, nchar(nm_curr)))
            fn_add <- ifelse(
              !is.null(value[[i]]$fn),
              paste0("_", value[[i]]$fn), "")
            out <- paste0(out, fn_add)
            params <- value[[i]]$params
            if (!is.null(params)) {
              params_add <- purrr::map_chr(seq_along(params), function(j) {
                nm_curr <- names(params)[j]
                val <- params[[j]]
                switch(
                  nm_curr,
                  "knots" = paste0(
                    "_k",
                    switch(
                      as.character(length(unique(diff(val)))),
                      "0" = "0",
                      "1" = paste0("e_", length(val)),
                      paste0("en_",  length(val))
                    )
                  ),
                  "Boundary.knots" = paste0("_b", length(val)),
                  "df" =  paste0("df", val)
                )
              })
              params_add <- paste0(params_add, collapse = "")
              if(length(params_add) > 0) out <- paste0(out, params_add)
            }
            out}),
            collapse = "_")
        )
      ),
    stop(paste0("name ", nm, " not recognised in calc_dir_proj"))
  )
  if(is.na(out) && !is.na(value)) {
    stop(paste0("NA returned for ", nm, " when value is ",
                paste0(value, sep = "/", collapse = "/")))
  }
  if(length(out) > 1) stop(paste0("Folder of length > 1 returned for ", nm,
                                  " when value is ",
                                  paste0(value, sep = "/", collapse = "/")))
  out
}

