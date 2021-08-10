# =========================
# data_raw
# =========================

#' @title Prepare raw data for the run
#' @export
prep_data_raw <- function(rmd, iter, ...) {

  .prep_data_raw <- switch(
    rmd,
    "cytokines" = .prep_dr_cytokines,
    "inf_markers" = .prep_dr_inf_markers,
    "hladr" = .prep_dr_hladr,
    "flowsom" = .prep_dr_flowsom,
    stop(paste0(rmd, " not recognised in prep_data_raw")))

  .prep_data_raw(iter = iter, ...)
}

.prep_dr_flowsom <- function(iter, stage, data_raw, ...) {

  switch(
    stage,
    "outer" = {

      # dataset
      data_raw <- switch(
        iter$ds,
        "cd4_th1_il17" = DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$flowsom,
        "cd8_th1" = DataTidyACSCyTOFCytokinesTCells::cd8_th1$flowsom,
        "tcrgd_th1" = DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$flowsom
      )

      # response type
      exc <- switch(
        iter$ds,
        "cd4_th1_il17" = "exc-Nd146Di",
        "cd8_th1" = "exc-Nd146Di",
        "tcrgd_th1" = "exc-Nd146Di"
      )

      # remove all_u- from stim
      data_raw <- data_raw %>%
        dplyr::mutate(stim = gsub("all_u-", "", stim))

      # filter
      filter_tbl <- iter[,c("stim", "non_responder_inclusion",
                           "n_clusters")]

      filter_tbl <- filter_tbl %>%
        dplyr::rename(
          n_clust = n_clusters,
          responders_only = non_responder_inclusion
        ) %>%
        dplyr::mutate(
          responders_only = ifelse(responders_only == "inc",
                                   "all", "r_o"),
          exc = exc)

      data_raw <- data_raw %>%
        filter_using_list(
          filter_list = filter_tbl
        )

      # remove clusters that have low stability
      data_raw <- data_raw %>%
        dplyr::filter(stability > 0.5)

      # remove non-responding subsets
      data_raw <- data_raw %>%
        dplyr::group_by(clust) %>%
        dplyr::select(SubjectID, VisitType, stim, exc, gn, scale, n_clust, clust,
                      n_cell_tot_stim, n_cell_tot_uns,
                      count_stim, count_uns, prob) %>%
        dplyr::filter(median(prob, na.rm = TRUE) > 0.5 & (sum(!is.na(prob)) / dplyr::n()) > 0.75) %>%
        dplyr::ungroup()

      # add in total number of cells in sample
      # for samples that are missing it
      # (missing presumably because no cyt+ cells were gated
      # out for these samples)
      n_cell_tbl <- switch(
        iter$ds,
        "cd4_th1_il17" = DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$stats_combn_tbl,
        "cd8_th1" = DataTidyACSCyTOFCytokinesTCells::cd8_th1$stats_combn_tbl,
        "tcrgd_th1" = DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$stats_combn_tbl
      )
      n_cell_tbl <- n_cell_tbl %>%
        dplyr::group_by(SubjectID, VisitType, stim) %>%
        dplyr::summarise(
          n_cell_tot_stim = n_cell_stim[1],
          n_cell_tot_uns = n_cell_uns[1],
          .groups = "drop"
        )

      data_raw <- data_raw %>%
        dplyr::select(-c(n_cell_tot_stim, n_cell_tot_uns)) %>%
        dplyr::left_join(
          n_cell_tbl,
          by = c("SubjectID", "VisitType", "stim")
        )

      # replace count_uns & count_stim with zero
      data_raw <- data_raw %>%
        dplyr::mutate(
          count_stim = ifelse(is.na(count_stim), 0, count_stim),
          count_unst = ifelse(is.na(count_uns), 0, count_uns),
        )

      # calculating freq_stim and freq_uns,
      # as in the per-cluster freq_stim and freq_stim
      data_raw <- data_raw %>%
        cytoutils::calc_freq(
          num = "count_stim",
          den = "n_cell_tot_stim",
          nm = "freq_tot_clust_stim"
        ) %>%
        cytoutils::calc_freq(
          num = "count_uns",
          den = "n_cell_tot_uns",
          nm = "freq_tot_clust_uns"
        )
      # calculate freq_bs, and then calculate the
      # number of ag + cells
      data_raw <- data_raw %>%
        dplyr::mutate(
          freq_tot_clust_bs = freq_tot_clust_stim - freq_tot_clust_uns,
          count_clust_ag = (freq_tot_clust_bs / 1e2) * n_cell_tot_stim
        )

      # calculate the number of ag-specific cells
      data_raw <- data_raw %>%
        dplyr::group_by(SubjectID, VisitType) %>%
        dplyr::mutate(count_clust_ag = ifelse(prob < 0.33, 0, count_clust_ag)) %>%
        dplyr::mutate(n_cell_ag = sum(count_clust_ag)) %>%
        dplyr::ungroup()

      # calculate the frequency of ag-specific
      # cells that belong to a given cluster
      data_raw <- data_raw %>%
        dplyr::mutate(freq_ag = count_clust_ag /
                        n_cell_ag * 1e2) %>%
        dplyr::mutate(freq_ag = pmax(freq_ag, 0)) %>%
        dplyr::mutate(count_ag = (freq_ag / 1e2) * n_cell_ag)

      # check that there is only one
      # observation per sample and cluster
      check_g_1 <- data_raw %>%
        dplyr::group_by(SubjectID, VisitType, clust) %>%
        dplyr::summarise(cnt = dplyr::n(),
                         .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::summarise(n_g_1 = any(cnt > 1)) %>%
        dplyr::pull(n_g_1)
      if(check_g_1) stop(
        "at least one key has more than one entry"
      )

      data_raw <- switch(
        iter$response_type,
        "freq_tot_bs" = data_raw %>%
          dplyr::group_by(SubjectID, VisitType,
                          clust,
                          n_cell_tot_stim,
                          n_cell_tot_uns) %>%
          dplyr::summarise(
            prob = mean(prob),
            n_cell = n_cell_tot_stim[1],
            resp = sum(count_clust_ag)
            ) %>%
          dplyr::ungroup(),
        "freq_ag" = data_raw %>%
          dplyr::rename(
            n_cell = n_cell_ag,
            resp = count_clust_ag,
            prob = prob
          ) %>%
          dplyr::mutate(n_cell = pmax(round(n_cell), 1)) %>%
          dplyr::mutate(resp = round(resp) / n_cell) %>%
          dplyr::group_by(SubjectID, VisitType) %>%
          dplyr::filter((1-prod(1-prob)) > 0.5) %>%
          dplyr::ungroup() %>%
          dplyr::select(SubjectID, VisitType, clust,
                        n_cell_tot_stim, prob, count_stim,
                        count_uns, n_cell, resp) %>%
          dplyr::filter(length(unique(clust)) > 1),
        "prob" = data_raw %>%
          dplyr::rename(
            resp = prob
          ) %>%
          dplyr::select(SubjectID, VisitType, clust, resp) %>%
          dplyr::group_by(clust) %>%
          dplyr::mutate(sd = sd(resp, na.rm = TRUE),
                           iqr = quantile(.data[["resp"]], na.rm = TRUE, 0.75) -
                             quantile(.data[["resp"]], na.rm = TRUE, 0.25)) %>%
          dplyr::filter(sd(resp, na.rm = TRUE) > 0.2,
                        iqr > 0.3) %>%
          dplyr::ungroup() %>%
          dplyr::select(-c(sd, iqr)),
        stop(paste0(iter$response_type, " not recognised"))
      )

    },
    "inner" = switch(
      iter$response_type,
      "freq_tot_bs" = ,
      "freq_ag" = ,
      "prob" = data_raw %>%
        dplyr::filter(.data$clust == iter$clust),
      stop(paste0(iter$response_type, " not recognised"))
    )
  )
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

.prep_dr_cytokines <- function(iter, stage, data_raw, ...) {

  switch(
    as.character(!grepl("prob$|fs$", iter$cyt_response_type)),
    "TRUE" = .prep_dr_cytokines_freq(
      iter = iter,
      stage = stage,
      data_raw = data_raw
    ),
    "FALSE" = switch(
      as.character(grepl("fs$", iter$cyt_response_type)),
      "FALSE" = .prep_dr_cytokines_prob(
        iter = iter,
        stage = stage,
        data_raw = data_raw)
      ,
      "TRUE" = .prep_dr_cytokines_scores(
        iter = iter,
        stage = stage,
        data_raw = data_raw
        )
    )
  )
}

.prep_dr_cytokines_freq <- function(iter,
                                    stage,
                                    data_raw) {
  switch(
    stage,
    "outer" = {
      data_raw <- switch(
        iter$ds,
        "cd4_th1_il17" =
          DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$stats_combn_tbl,
        "cd8_th1" =
          DataTidyACSCyTOFCytokinesTCells::cd8_th1$stats_combn_tbl,
        "tcrgd_th1" =
          DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$stats_combn_tbl
      )

      # select required columns
      data_raw <- data_raw %>%
        dplyr::select(SubjectID, VisitType, stim, cyt_combn,
                      count_stim, n_cell_stim, count_uns,
                      n_cell_uns)

      # filter
      data_raw <- data_raw %>%
        filter_using_list(
          filter_list = iter["stim"]
        )

      # calculate background-subtracted frequencies, if need be
      data_raw <- data_raw %>%
        .subtract_background()


      # convert cyt combn from COMPASS format
      data_raw$cyt_combn <- compassutils::convert_cyt_combn_format(
        cyt_combn = data_raw$cyt_combn,
        to = "std"
      )

      # remove all-neg cyt_combn
      data_raw <- data_raw %>%
        dplyr::filter(stringr::str_detect(cyt_combn, "\\+"))

      # sum over markers
      data_raw <- switch(
        iter$cyt_response_type,
        "summed" = data_raw %>%
          cytoutils::sum_over_markers(
            grp = c("SubjectID", "VisitType",
                    "stim", "n_cell"),
            cmbn = "cyt_combn",
            levels = c("-", "+"),
            resp = c("freq_bs")
          ) %>%
          dplyr::mutate(cyt_combn = "summed") %>%
          dplyr::mutate(resp = freq_bs * n_cell),
        "cyt_prop" = ,
        "cyt" = {
          marker_vec <- stringr::str_split(unique(data_raw$cyt_combn),
                                           pattern = "\\+|\\-")[[1]]
          marker_vec <- marker_vec[-which(marker_vec == "")]
          purrr::map_df(marker_vec, function(mk) {
            data_raw %>%
              cytoutils::sum_over_markers(
                grp = c("SubjectID", "VisitType",
                        "stim", "n_cell"),
                cmbn = "cyt_combn",
                markers_to_sum = setdiff(marker_vec, mk),
                levels = c("-", "+"),
                resp = c("freq_bs", "resp")
              ) %>%
              dplyr::filter(grepl("\\+", cyt_combn))
          })
        },
        "combn" = ,
        "combn_prop" = data_raw,
        stop(paste0(iter$cyt_response_type, " not"))
      )

      # remove non-responding cyt_combn's
      # from analysis if not summed
      data_raw <- switch(
        iter$cyt_response_type,
        "summed" = data_raw,
        "cyt_prop" = ,
        "cyt" = {
          compass_obj <- switch(
            iter$ds,
            "cd4_th1_il17" =
              DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$compass,
            "cd8_th1" =
              DataTidyACSCyTOFCytokinesTCells::cd8_th1$compass,
            "tcrgd_th1" =
              DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$compass
          )  %>%
            magrittr::extract2("locb0.15_min_clust") %>%
            magrittr::extract2(iter$stim)

          if(is.null(compass_obj)) return(data_raw[1,][-1,])

          compass_obj <- compass_obj %>%
            magrittr::extract2("fit")%>%
            magrittr::extract2("mean_gamma")

          .remove_cyt_low_prob(
            .data = data_raw,
            compass_obj = compass_obj,
            quant_min = 0.25,
            prob_min = 0.8
          )
        },
        "combn_prop" = ,
        "combn" = {
          compass_obj <- switch(
            iter$ds,
            "cd4_th1_il17" =
              DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$compass,
            "cd8_th1" =
              DataTidyACSCyTOFCytokinesTCells::cd8_th1$compass,
            "tcrgd_th1" =
              DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$compass
          ) %>%
            magrittr::extract2("locb0.15_min_clust") %>%
            magrittr::extract2(iter$stim)

          if(is.null(compass_obj)) return(data_raw[1,][-1,])

          compass_obj <- compass_obj %>%
            magrittr::extract2("fit")%>%
            magrittr::extract2("mean_gamma")

          .remove_combn_low_prob(
            .data = data_raw,
            compass_obj = compass_obj,
            quant_min = 0.25,
            prob_min = 0.8
          )
        },
        stop(paste0(iter$cyt_response_type, " not recognised"))
      )

      if(nrow(data_raw) == 0) return(data_raw)

      # remove individuals that didn't respond to anything
      data_raw <- switch(
        iter$cyt_response_type,
        "summed" = ,
        "cyt" = ,
        "combn" = data_raw,
        "cyt_prop" = ,
        "combn_prop" = {
          prob_tbl <- switch(
            iter$ds,
            "cd4_th1_il17" =
              DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$post_probs_bulk %>%
              magrittr::extract2("exc-Nd146Di"),
            "cd8_th1" =
              DataTidyACSCyTOFCytokinesTCells::cd8_th1$post_probs_bulk %>%
              magrittr::extract2("exc-Nd146Di"),
            "tcrgd_th1" =
              DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$post_probs_bulk %>%
              magrittr::extract2("exc-Nd146Di")
          ) %>%
            magrittr::extract2("locb0.15_min_clust") %>%
            magrittr::extract2(iter$stim)

          if(is.null(prob_tbl)) return(data_raw[1,][-1,])

          data_raw %>%
            dplyr::filter(
              paste0(SubjectID, "_", VisitType) %in%
                prob_tbl[["sampleid"]][prob_tbl[["prob"]] > 0.75]
            )
        }
      )


      # sum or calculate proportion (or do nothing)
      data_raw <- switch(
        iter$cyt_response_type,
        "summed" = ,
        "cyt" = ,
        "combn" = data_raw,
        "cyt_prop" = ,
        "combn_prop" = {
          data_raw %>%
            dplyr::group_by(SubjectID, VisitType, stim) %>%
            dplyr::mutate(n_cell_ag = sum(resp)) %>%
            dplyr::mutate(n_cell_ag = round(n_cell_ag)) %>%
            dplyr::filter(n_cell_ag >= 1) %>%
            dplyr::mutate(resp = resp/n_cell_ag) %>%
            dplyr::mutate(resp = pmin(1, resp) %>%
                            pmax(0)) %>%
            dplyr::select(-c(n_cell, freq_bs)) %>%
            dplyr::rename(n_cell = n_cell_ag) %>%
            dplyr::ungroup() %>%
            dplyr::select(SubjectID, VisitType,
                          stim, cyt_combn, n_cell, resp) %>%
            dplyr::filter(length(unique(cyt_combn)) > 1) %>%
            dplyr::select(SubjectID, VisitType, cyt_combn, n_cell, resp) %>%
            dplyr::group_by(cyt_combn) %>%
            dplyr::mutate(sd = sd(resp, na.rm = TRUE),
                          iqr = quantile(.data[["resp"]], na.rm = TRUE, 0.75) -
                            quantile(.data[["resp"]], na.rm = TRUE, 0.25)) %>%
            dplyr::filter(sd(resp, na.rm = TRUE) > 0.125,
                          iqr > 0.2) %>%
            dplyr::ungroup() %>%
            dplyr::select(-c(sd, iqr))
        },
        stop(paste0(iter$cyt_response_type, " not recognised"))
      )
    },
    "inner" = switch(
      iter$cyt_response_type,
      "summed" = data_raw,
      "cyt" = ,
      "cyt_prob" = ,
      "cyt_prop" = ,
      "combn" = ,
      "combn_prob" = ,
      "combn_prop" = {
        combn_comp <- paste0(
          "^",
          paste0(cytoutils:::add_double_backslash(
            iter$cyt_response_type_grp
          )),
          "$")
        data_raw <- data_raw %>%
          dplyr::filter(grepl(combn_comp,
                              .data$cyt_combn))

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

      }
    )
  )
}

.prep_dr_cytokines_prob <- function(iter,
                                    stage,
                                    data_raw) {

  switch(
    stage,
    "outer" = {
      data_raw <- switch(
        iter$ds,
        "cd4_th1_il17" =
          DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$compass,
        "cd8_th1" =
          DataTidyACSCyTOFCytokinesTCells::cd8_th1$compass,
        "tcrgd_th1" =
          DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$compass
      ) %>%
        magrittr::extract2("locb0.15_min_clust") %>%
        magrittr::extract2(iter$stim)
      if(is.null(data_raw)) return(
        tibble::tibble(SubjectID = character(0),
                       VisitType = character(0),
                       cyt_combn = character(0),
                       resp = numeric(0))
      )

      data_raw <- data_raw %>%
        magrittr::extract2("fit")%>%
        magrittr::extract2("mean_gamma")

      rn_vec <- rownames(data_raw)

      data_raw <- data_raw %>%
        tibble::as_tibble() %>%
        dplyr::mutate(SampleID = rn_vec) %>%
        dplyr::select(SampleID, everything())

      colnames(data_raw)[-1] <-
        compassutils::convert_cyt_combn_format(
          colnames(data_raw)[-1], to = "std"
        )

      data_raw <- data_raw %>%
        tidyr::pivot_longer(
          -SampleID,
          names_to = "combn",
          values_to = "prob"
        )

      # remove all-neg cyt_combn
      data_raw <- data_raw %>%
        dplyr::filter(stringr::str_detect(combn, "\\+"))

      # calculate per-cyt/cyt_combn probability of a response
      data_raw <- switch(
        iter$cyt_response_type,
        "combn_prob" = data_raw,
        "cyt_prob" = {
          marker_vec <- stringr::str_split(unique(data_raw$combn),
                                           pattern = "\\+|\\-")[[1]]
          marker_vec <- marker_vec[-which(marker_vec == "")]
          purrr::map_df(marker_vec, function(mk) {
            data_raw %>%
              dplyr::filter(grepl(
                paste0(mk, "\\+"),
                combn
              )) %>%
              dplyr::mutate(combn = paste0(mk, "+")) %>%
              dplyr::group_by(SampleID, combn) %>%
              dplyr::summarise(prob = 1 - prod(1 - prob),
                               .groups = "drop_last") %>%
              dplyr::ungroup()
          })

        }
      )

      quant_min <- 0.25; prob_min <- 0.8
      combn_vec_sel <- data_raw %>%
        dplyr::group_by(combn) %>%
        dplyr::filter(quantile(.data$prob, 1 - quant_min) >= prob_min) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::pull("combn")

      data_raw <- data_raw %>%
        dplyr::filter(combn %in% combn_vec_sel)

      data_raw %>%
        tidyr::separate(
          col = SampleID,
          into = c("SubjectID", "VisitType"),
          sep = "_"
        ) %>%
        dplyr::rename(cyt_combn = combn) %>%
        dplyr::rename(resp = prob)
    },
    "inner" = switch(
      iter$cyt_response_type,
      "cyt_prob" = ,
      "combn_prob" = {
        combn_comp <- paste0(
          "^",
          paste0(cytoutils:::add_double_backslash(
            iter$cyt_response_type_grp
          )),
          "$")
        data_raw <- data_raw %>%
          dplyr::filter(grepl(combn_comp,
                              .data$cyt_combn))

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

      }
    )
  )

}
.prep_dr_cytokines_scores <- function(iter,
                                      stage,
                                      data_raw) {

  switch(
    stage,
    "outer" = {
      data_raw <- switch(
        iter$ds,
        "cd4_th1_il17" =
          DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$compass,
        "cd8_th1" =
          DataTidyACSCyTOFCytokinesTCells::cd8_th1$compass,
        "tcrgd_th1" =
          DataTidyACSCyTOFCytokinesTCells::tcrgd_th1$compass
      ) %>%
        magrittr::extract2("locb0.15_min_clust") %>%
        magrittr::extract2(iter$stim)

      if(is.null(data_raw)) return(
        tibble::tibble(SubjectID = character(0),
                       VisitType = character(0),
                       cyt_combn = character(0),
                       resp = numeric(0))
      )

      data_raw <- data_raw %>%
        COMPASS::scores() %>%
        tibble::as_tibble() %>%
        dplyr::select(SampleID, SubjectID, VisitType, FS, PFS) %>%
        dplyr::rename(fs = FS,
                      pfs = PFS)

      data_raw <- data_raw[,-which(colnames(data_raw) ==
                                     setdiff(c("fs", "pfs"),
                                             iter$cyt_response_type))]

      colnames(data_raw)[which(colnames(data_raw) ==
                        iter$cyt_response_type)] <- "resp"

      data_raw %>%
        dplyr::mutate(cyt_combn = iter$cyt_response_type)

    },
    "inner" = switch(
      iter$cyt_response_type,
      "fs" = ,
      "pfs" = data_raw
    )
  )
}

.remove_combn_low_prob <- function(.data, compass_obj, quant_min, prob_min) {

  compass_mat_prob <- compass_obj
  rn_vec <- rownames(compass_mat_prob)

  compass_tbl_prob <- compass_mat_prob %>%
    tibble::as_tibble() %>%
    dplyr::mutate(SampleID = rn_vec) %>%
    dplyr::select(SampleID, everything())

  colnames(compass_tbl_prob)[-1] <-
    compassutils::convert_cyt_combn_format(
      colnames(compass_tbl_prob)[-1], to = "std"
    )

  compass_tbl_prob <- compass_tbl_prob %>%
    magrittr::extract(,c("SampleID", unique(.data$cyt_combn)))

  combn_vec_sel <- compass_tbl_prob %>%
    tidyr::pivot_longer(-SampleID,
                        names_to = "combn",
                        values_to = "prob") %>%
    dplyr::group_by(combn) %>%
    dplyr::filter(quantile(.data$prob, 1 - quant_min) >= prob_min) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull("combn")

  .data %>%
    dplyr::filter(cyt_combn %in% combn_vec_sel)

}

.remove_cyt_low_prob <- function(.data, compass_obj, quant_min, prob_min) {

  compass_mat_prob <- compass_obj
  rn_vec <- rownames(compass_mat_prob)

  compass_tbl_prob <- compass_mat_prob %>%
    tibble::as_tibble() %>%
    dplyr::mutate(SampleID = rn_vec) %>%
    dplyr::select(SampleID, everything())

  colnames(compass_tbl_prob)[-1] <-
    compassutils::convert_cyt_combn_format(
      colnames(compass_tbl_prob)[-1], to = "std"
    )

  compass_tbl_prob <- compass_tbl_prob %>%
    tidyr::pivot_longer(
      -SampleID,
      names_to = "cyt_combn",
      values_to = "prob"
    )

  marker_vec <- unique(.data$cyt_combn)
  marker_vec <- stringr::str_remove_all(marker_vec, "\\+|\\-")
  compass_tbl_prob <- purrr::map_df(marker_vec, function(mk) {
    compass_tbl_prob %>%
      dplyr::filter(grepl(
        paste0(mk, "\\+"),
        cyt_combn
      )) %>%
      dplyr::mutate(combn = paste0(mk, "+")) %>%
      dplyr::group_by(SampleID, combn) %>%
      dplyr::summarise(prob = 1 - prod(1 - prob),
                       .groups = "drop_last") %>%
      dplyr::ungroup()
  })

  combn_vec_sel <- compass_tbl_prob %>%
    dplyr::group_by(combn) %>%
    dplyr::filter(quantile(.data$prob, 1 - quant_min) >= prob_min) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull("combn")

  .data %>%
    dplyr::filter(cyt_combn %in% combn_vec_sel)


}

.subtract_background <- function(.data) {
  .data %>%
    cytoutils::calc_prop(
      den = "n_cell_stim",
      num = "count_stim",
      nm = "prop_stim"
    ) %>%
    cytoutils::calc_prop(
      den = "n_cell_uns",
      num = "count_uns",
      nm = "prop_uns"
    ) %>%
    dplyr::mutate(
      prop_bs = prop_stim - prop_uns,
      count_bs = prop_bs * n_cell_stim
    ) %>%
    dplyr::rename(
      n_cell = n_cell_stim,
      resp = count_bs
    ) %>%
    dplyr::rename(freq_bs = prop_bs) %>%
    dplyr::select(-c(count_stim, count_uns,
                     n_cell_uns, prop_stim,
                     prop_uns))
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
    "cytokines" = .prep_iter_cytokines,
    "hladr" = .prep_iter_identity,
    "inf_markers" = .prep_iter_identity,
    "flowsom" = .prep_iter_flowsom,
    stop(paste0("no .prep_iter fn defined for ", rmd))
  )
  .prep_iter(iter = iter, ...)
}

.prep_iter_identity <- function(iter, ...) iter

.prep_iter_cytokines <- function(iter, data_raw) {

  purrr::map_df(unique(data_raw$cyt_combn),
                function(cyt_combn){
                  iter %>%
                    dplyr::mutate(
                      cyt_response_type_grp = cyt_combn)
                })
}

.prep_iter_flowsom <- function(iter, data_raw) {

  if(nrow(data_raw) == 0) return(iter)

  purrr::map_df(unique(data_raw$clust),
                function(clust){
                  iter %>%
                    dplyr::mutate(
                      clust = clust
                      )
                })
}

