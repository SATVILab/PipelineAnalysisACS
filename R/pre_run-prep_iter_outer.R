# =========================
# prep iter_outer
# =========================

#' @title Prepare outer iterator
#'
#' @description Add or filter it, as needed.
#'
#' @export
prep_iter <- function(rmd, iter, ...) {
  if (is.null(names(iter$var_exp_spline))) {
    iter$var_exp_spline <- iter$var_exp_spline[[1]]
  }
  .prep_iter <- switch(rmd,
    "cytokines" = .prep_iter_cytokines,
    "hladr" = .prep_iter_identity,
    "inf_markers" = .prep_iter_identity,
    "flowsom" = .prep_iter_flowsom,
    "faust" = .prep_iter_faust,
    "faust_cyt" = .prep_iter_faust_cyt,
    stop(paste0("no .prep_iter fn defined for ", rmd))
  )
  .prep_iter(iter = iter, ...)
}

.prep_iter_faust <- function(iter, data_raw) {
  purrr::map_df(unique(data_raw$den), function(den) {
    data_raw_den <- data_raw %>%
      dplyr::filter(den == .env$den)
    purrr::map_df(unique(data_raw_den$pop_sub_faust), function(pop_sub_faust) {
      iter %>%
        dplyr::mutate(den = den, pop_sub_faust = pop_sub_faust)
    })
  })
}

.prep_iter_faust_cyt <- function(iter, data_raw, ...) {
  dr_add_col_vec <- setdiff(
    c("pop", "pheno", "combn", "stim"),
    colnames(iter)
  )

  if (nrow(data_raw) == 0) {
    return(iter)
  }

  data_raw_add <- data_raw %>%
    dplyr::select_at(dr_add_col_vec) %>%
    dplyr::group_by_at(dr_add_col_vec) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  purrr::map_df(seq_len(nrow(iter)), function(i) {
    iter_bind <- purrr::map_df(seq_len(nrow(data_raw_add)), function(x) {
      iter[i, ]
    })
    iter_bind %>%
      dplyr::bind_cols(
        data_raw_add
      )
  })
}

.prep_iter_identity <- function(iter, ...) iter

.prep_iter_cytokines <- function(iter, data_raw, ...) {
  purrr::map_df(
    unique(data_raw$cyt_combn),
    function(cyt_combn) {
      iter %>%
        dplyr::mutate(
          cyt_response_type_grp = cyt_combn
        )
    }
  )
}

.prep_iter_flowsom <- function(iter, data_raw, ...) {
  if (nrow(data_raw) == 0) {
    return(iter)
  }

  purrr::map_df(
    unique(data_raw$clust),
    function(clust) {
      iter %>%
        dplyr::mutate(
          clust = clust
        )
    }
  )
}
