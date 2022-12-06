#' @title Get project directory based on iterators
#'
#' @export
calc_dir_proj <- function(dir_base, iter = list(), ...) {
  d_list <- iter |> append(list(...))
  dir_sub <- purrr::map_chr(
    seq_along(d_list),
    function(i) {
      out <- get_folder(
        nm = names(d_list)[i],
        value = d_list[[i]]
      )
      out
    }
  ) |>
    paste0(collapse = "/")
  file.path(dir_base, dir_sub) |>
    fs::path_norm()
}

get_folder <- function(nm, value) {
  value <- switch(as.character(is.na(value)),
    "TRUE" = "-",
    "FALSE" = value
  )
  out <- switch(nm,
    "ttb_max" = as.character(value),
    "var_conf" = ,
    "conf" = ifelse(
      identical(value[[1]], "none"), "no_conf", paste0(
        stringr::str_sub(value[[1]], end = min(nchar(value[[1]]), 3)),
        collapse = ""
      )
    ),
    "filter_approach" = switch(names(value[[1]]),
      "filter" = {
        out <- "fltr/"
        f_l <- value[[1]][[1]]
        for (i in seq_along(f_l)) {
          nm_elem <- names(f_l)[i]
          val_elem <- f_l[[i]]
          val_elem <- paste0(names(val_elem), val_elem, collapse = "_")
          out <- paste0(
            nm_elem, "_", val_elem
          )
        }
        out
        # override above
        "filter"
      },
      "select" = {
        outer <- "slct/"
        # override above
        "select"
      },
      stop("names(iter$filter_approach) not recognised")
    ),
    "den" = value,
    "pop" = value,
    "combn" = {
      out <- ""
      cyt_vec <- c("IFNg-beads", "IL2", "IL6", "IL17", "IL22", "TNFa")
      for (cyt in cyt_vec) {
        pos <- stringr::str_locate(value, cyt)
        ind <- stringr::str_sub(
          value,
          pos[, "end"] + 1,
          pos[, "end"] + 1
        )
        if (grepl("\\+", ind)) {
          out <- paste0(out, cyt, "+")
        }
      }
      out
    },
    "pheno" = {
      full_pheno_label <- "CD4~2~2~HLA-DR-beads~1~2~CD8-IgD~1~3~CD3~2~2~CD7~2~2~CD27~1~2~TCRgd-CD19~1~2~CD33~1~2~CD16~1~2~CCR7~1~2~CXCR5~1~2~Perf-beads~1~2~CD45RA~1~2~CD14~1~2~CD161~1~2~CD20~1~2~"
      marker_vec <- stringr::str_split(full_pheno_label, "~\\d~\\d~")[[1]]
      marker_vec <- marker_vec[marker_vec != ""]
      out <- ""
      for (marker in marker_vec) {
        if (!grepl(marker, value)) next
        pos <- stringr::str_locate(value, marker)
        ind <- stringr::str_sub(
          value,
          pos[, "end"] + 1,
          pos[, "end"] + 1
        )
        if (grepl("\\+", ind)) {
          out <- paste0(out, marker, "+")
        }
      }
      out <- ifelse(out == "", "all-", out)
      out
    },
    "pkg" = value,
    "rem_il6" = ifelse(value, "il6_e", "il6_i"),
    "wins" = switch(as.character(is.logical(value)),
      "TRUE" = ifelse(value, "wins_y", "wins_n"),
      "FALSE" = value
    ),
    "equi_d_knots" = ifelse(value, "eq_d", "eq_d_n"),
    "den" = value,
    "pop_sub_faust" = value,
    "ds" = ,
    "dataset" = value,
    "stim" = value,
    "var_dep" = value,
    "cyt_response_type" = value,
    "var_offset" = paste0(
      "o-", stringr::str_sub(value, end = min(nchar(value), 3)
      ) |>
      paste0(collapse = "_")),
    "cyt_response_type_grp" = value,
    "response_type" = value,
    "non_responder_inclusion" = paste0("n_r-", value),
    "n_clusters" = value,
    "clust" = value,
    "var_exp" = paste0("exp-", paste0(
      stringr::str_sub(value, end = min(nchar(value), 5)),
      collapse = ""
    )),
    "var_int" = ifelse(
      is.null(value[[1]]),
      "int_n",
      "int"
    ),
    "mod" = value,
    "var_exp_spline" =
      switch(as.character(identical(value, "none")),
        "TRUE" = "exp_s-none",
        "FALSE" = paste0(

          # "exp_s-",
          "",
          paste0(purrr::map_chr(seq_along(value), function(i) {
            nm_curr <- names(value)[i]
            if (grepl("^tc~", nm_curr)) {
              nm_curr <- stringr::str_split(nm_curr, "~")[[1]][3]
            }
            out <- stringr::str_sub(nm_curr, end = min(6, nchar(nm_curr)))
            fn_add <- ifelse(
              !is.null(value[[i]]$fn),
              paste0("_", value[[i]]$fn), ""
            )
            out <- paste0(out, fn_add)
            params <- value[[i]]$params
            if (!is.null(params)) {
              params_add <- purrr::map_chr(seq_along(params), function(j) {
                nm_curr <- names(params)[j]
                val <- params[[j]]
                switch(nm_curr,
                  "knots" = paste0(
                    "_k",
                    switch(as.character(length(unique(diff(val)))),
                      "0" = "0",
                      "1" = paste0("e_", length(val)),
                      paste0("en_", length(val))
                    )
                  ),
                  "Boundary.knots" = paste0("_b", length(val)),
                  "df" = paste0("df", val)
                )
              })
              params_add <- paste0(params_add, collapse = "")
              if (length(params_add) > 0) out <- paste0(out, params_add)
            }
            out
          }),
          collapse = "_"
          )
        )
      ),
    stop(paste0("name ", nm, " not recognised in calc_dir_proj"))
  )
  if (is.na(out) && !is.na(value)) {
    stop(paste0(
      "NA returned for ", nm, " when value is ",
      paste0(value, sep = "/", collapse = "/")
    ))
  }
  if (length(out) > 1) {
    stop(paste0(
      "Folder of length > 1 returned for ", nm,
      " when value is ",
      paste0(value, sep = "/", collapse = "/")
    ))
  }
  out
}
