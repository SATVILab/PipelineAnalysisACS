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
  .prep_iter <- switch(
    rmd,
    "cyt_pos_freq" = .prep_iter_cyt_pos_freq,
    "hladr" = .prep_iter_identity,
    "inf_markers" = .prep_iter_identity,
    stop(paste0("no .prep_iter fn defined for ", rmd))
  )
  .prep_iter(iter = iter, ...)
}

.prep_iter_identity <- function(iter, ...) iter

.prep_iter_cyt_pos_freq <- function(iter, data_raw) {

  purrr::map_df(unique(data_raw$cyt_combn),
                function(cyt_combn){
                  iter %>%
                    dplyr::mutate(
                      cyt_response_type_grp = cyt_combn)
                })
}

# =========================
# data_raw
# =========================

#' @title Prepare raw data for the run
#' @export
prep_data_raw <- function(rmd, iter, ...) {

  .prep_data_raw <- switch(
    rmd,
    "cyt_pos_freq" = .prep_dr_cyt_pos_freq,
    "inf_markers" = .prep_dr_inf_markers,
    "hladr" = .prep_dr_hladr,
    stop(paste0(rmd, " not recognised in prep_data_raw")))

  .prep_data_raw(iter = iter, ...)
}

.prep_dr_hladr <- function(iter, stage, data_raw, ...) {

  switch(
    stage,
    "outer" = {
      data_raw <- DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$hladr

      data_raw <- data_raw %>%
        filter_using_list(
          filter_list = iter["stim"]
        )

      data_raw <- data_raw %>%
        dplyr::mutate(
          SubjectID = substr(.data$fcs, start = 1, stop = 6),
          fcs = substr(.data$fcs, start = 8, stop = nchar(fcs)),
          VisitType = purrr::map_chr(fcs, function(x) {
            uns_loc <- gregexpr("_", x)[[1]][1]
            substr(x, start = 1, stop = uns_loc - 1)
          }))

      data_raw <- data_raw %>%
        dplyr::rename(resp = hladr_med_diff)

      data_raw <- data_raw %>%
        dplyr::select(-c(batch, batch_sh, fcs, ind, is_uns, mult,
                         ind_in_batch, hladr_med_pos, hladr_med_neg,
                         chnl, gate_name)) %>%
        dplyr::select(SubjectID, VisitType, stim, resp)

      check_g_1 <- data_raw %>%
        dplyr::group_by(SubjectID, VisitType) %>%
        dplyr::summarise(cnt = dplyr::n(),
                         .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::summarise(n_g_1 = any(cnt > 1)) %>%
        dplyr::pull(n_g_1)
      if(check_g_1) stop(
        "at least one key has more than one entry"
      )

      data_raw

    },
    "inner" = data_raw
  )
}

.prep_dr_cyt_pos_freq <- function(iter, stage, data_raw, ...) {

  switch(
    stage,
    "outer" = {
      data_raw <- switch(
        iter$ds,
        "cd4_th1_il17" = DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$stats_combn_tbl,
        "cd8_th1" = DataTidyACSCyTOFCytokinesTCells::cd8_th1$stats_combn_tbl,
        "tcrgd_th1" = DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$stats_combn_tbl
      )

      # filter
      data_raw <- data_raw %>%
        filter_using_list(
          filter_list = iter["stim"]
        )

      data_raw <- prep_bs_freq(
        .data = data_raw,
        cyt_response_type = iter$cyt_response_type
      )
    },
    "inner" = data_raw
  )

}

.prep_dr_inf_markers <- function(iter, stage, data_raw, ...) {

  switch(
    stage,
    "outer" =   switch(
      iter$var_dep,
      "risk6" = TuberculomicsCompendium::signature_6gene %>%
        tibble::as_tibble() %>%
        dplyr::rename(resp = sig6gene_CorScore),
      TuberculomicsCompendium::soma_data_tidy %>%
        dplyr::filter(Soma_Target == iter$var_dep) %>%
        dplyr::group_by_at(c("SubjectID", "VisitType", "Soma_Target")) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::rename(resp = Soma_TransformedReadout) %>%
        dplyr::select(SubjectID, VisitType, resp)
    ),
    "inner" = data_raw
  )

}

