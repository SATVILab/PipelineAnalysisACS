#' @export
filter_using_list <- function(.data, filter_list) {
  for (i in seq_along(filter_list)) {
    nm <- names(filter_list)[i]
    value <- filter_list[[i]]
    .data <- .data[.data[[nm]] %in% value,]
  }
  .data
}

trans <- function(.data, trans) {
  if(is.null(trans)) return(.data)
  if(is.function(trans)) {
    .data <- .data %>%
      dplyr::mutate(resp = trans(resp))
    return(.data)
  }
  if (!is.character(trans)) {
    stop("trans must be NULL, a function or a character")
  }

  trans_fn <- switch(
    trans,
    "within_01" = function(x) {
      range_x <- diff(range(x))
      pmax(pmin(max(x) - range_x * 0.01, x), range_x * 0.01)
    },
    "pos" = function(x) pmax(0, x),
    "pos_strict" = function(x) {
      range_x <- diff(range(x))
      pmax(range_x * 0.01, x)
    },
    "sqrt" = sqrt,
    "sqrt_pos" = function(x) sqrt(pmax(0, x)),
    stop("trans supplied but not recognised")
  )
  .data %>%
    dplyr::mutate(resp = trans_fn(resp))
}

#' @title Calculate bs-freq for CyTOF phenotypes
#'
#' @param .data dataframe.
#' @param cyt_response_type "summed" or "combn". Specifies
#' type of response to create.
#' Should be able to specify "single" as well.
#' @export
prep_bs_freq <- function(.data, cyt_response_type) {

  .data$cyt_combn <- compassutils::convert_cyt_combn_format(
    cyt_combn = .data$cyt_combn,
    to = "std"
  )
  .data <- .data %>%
    dplyr::filter(grepl("\\+", .data$cyt_combn))

  # sum over cytokines
  .data <- switch(
    cyt_response_type,
    "summed" = .data %>%
      cytoutils::sum_over_markers(
        grp = c("batch_sh", "SubjectID", "VisitType",
                "stim", "n_cell_stim", "n_cell_uns"),
        cmbn = "cyt_combn",
        levels = c("-", "+"),
        resp = c("count_stim", "count_uns")
      ) %>%
      dplyr::mutate(cyt_combn = "summed"),
    "combn" = .data %>%
      dplyr::select(batch_sh, SubjectID, VisitType,
                    cyt_combn,
                    stim, n_cell_stim, n_cell_uns,
                    count_stim, count_uns),
    stop(paste0(
        paste0(cyt_response_type, collapse = "/"),
        " not available for cyt_response_type in prep_bs_freq")
        )
    )

  .data <- .data %>%
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
    dplyr::mutate(prop_bs = pmax(prop_bs * 1e2, 0),
                  resp = pmax(resp, 0)) %>%
    dplyr::rename(freq_bs = prop_bs) %>%
    dplyr::select(-c(count_stim, count_uns,
                     n_cell_uns, prop_stim,
                     prop_uns))

  .data
}

#' @export
#'
#' @param
#'
#' @details
#' - Join on SubjectID and VisitType
#' - Always add Progressor, tfmttb and timeToTBFromVisit
#' - Returns only rows from .data
#' where sample is before TB diagnosis
add_clinical_data_and_filter <- function(.data,
                              cols_add = NULL,
                              cols_join = c("SubjectID", "VisitType"),
                              cols_add_perm =  c("Progressor",
                                                 "timeToTBFromVisit",
                                                 "tfmttb"),
                              ttb_min = NULL,
                              ttb_max) {

  cols_add <- unique(c(cols_add_perm, cols_add)) %>%
    setdiff("n_cell")
  cols_add <- unique(cols_add)
  cols_join <- unique(cols_join)
  clinical_data_add <- TuberculomicsCompendium::clinical_data
  clinical_data_add <- clinical_data_add[,colnames(clinical_data_add) %in%
                                           c(cols_join, cols_add)]
  clinical_data_add <- clinical_data_add %>%
    dplyr::mutate(SampleID = paste0(.data$SubjectID, "_", .data$VisitType)) %>%
    dplyr::filter(SampleID %in% unique(
      .env$.data %>%
        dplyr::mutate(SampleID = paste0(.data$SubjectID, "_", .data$VisitType)) %>%
        dplyr::pull("SampleID")
    )) %>%
    dplyr::select(-SampleID)

  clinical_data_add <- clinical_data_add %>%
    dplyr::group_by_at(c(cols_join, setdiff(cols_add, "tfmttb"))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # check that there is only one row per key
  # (key is combn of joining col levels)
  clinical_data_add_check <- clinical_data_add %>%
    dplyr::group_by_at(cols_join) %>%
    dplyr::summarise(cnt = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(g_1 = any(cnt > 1)) %>%
    dplyr::pull(g_1)
  if(clinical_data_add_check) stop("clinical_data has more than one row per key")

  clinical_data_add <- clinical_data_add %>%
    dplyr::filter(timeToTBFromVisit >= 0 | is.na(timeToTBFromVisit))

  clinical_data_add <- clinical_data_add %>%
    dplyr::filter(Progressor %in% c("yes", "no"))

  if (!is.null(ttb_min)) {
    clinical_data_add <- clinical_data_add %>%
      dplyr::mutate(timeToTBFromVisit = pmax(.data$timeToTBFromVisit, ttb_min))
  }

  clinical_data_add <- clinical_data_add %>%
    dplyr::mutate(tfmttb = ifelse(Progressor == "no", # make zero non-progressor
                                  0,
                                  -9999),
                  tfmttb = ifelse(.data$Progressor == "yes", # make 0 if time to TB is greater than ttb_max; otherwise make
                                  # how many days it is from collapse time to the actual time
                                  pmax(0, .env$ttb_max - .data$timeToTBFromVisit),
                                  .data$tfmttb))

  .data <- .data[,!(colnames(.data) %in%
                      setdiff(colnames(clinical_data_add),
                              cols_join))]

  cn_vec_orig <- colnames(.data) %>%
    setdiff(cols_join)

  .data <- .data %>%
    dplyr::inner_join(
      clinical_data_add,
      by = cols_join
    )

 .data <- .data[,c(cols_join,
                   cols_add,
                   setdiff(cn_vec_orig,
                           c(cols_join, cols_add))) %>%
                  unique()]

 .data

}

#' @title Add assay data
add_tc_assay_data <- function(.data,
                              assay,
                              cols) {


  if(any(purrr::map_lgl(p_dots$conf, function(x) stringr::str_detect(x, "sig6gene_CorScore")))){
    data_raw %<>%
      dplyr::left_join(TuberculomicsCompendium::signature_6gene,
                       by = c("SubjectID", "VisitType")) %>%
      dplyr::select(dataset:Winter, sig6gene_CorScore, PreviousDiagnosisOfTB,
                    everything())
    data_raw %<>%
      dplyr::filter(!is.na(sig6gene_CorScore))
  }

  if(names(p_dots$var_exp_spline) != "tfmttb"){
    if(names(p_dots$var_exp_spline) %in% TuberculomicsCompendium::soma_data_tidy$Soma_Target){
      var_tbl_add <- TuberculomicsCompendium::soma_data_tidy %>%
        dplyr::filter(Soma_Target == names(p_dots$var_exp_spline)) %>%
        dplyr::group_by(SubjectID, VisitType) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(SubjectID, VisitType, Soma_TransformedReadout)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(p_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add ,
                   by = c("SubjectID", "VisitType"))
    } else if(names(p_dots$var_exp_spline) == "risk6"){
      var_tbl_add <- TuberculomicsCompendium::signature_6gene %>%
        dplyr::select(SubjectID, VisitType, sig6gene_CorScore)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(p_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType"))
    } else if(names(p_dots$var_exp_spline) == "risk11"){
      var_tbl_add <- TuberculomicsCompendium::signature_11gene %>%
        dplyr::select(SubjectID, VisitType, sig11gene_CorScore)
      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(p_dots$var_exp_spline))
      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType"))
    } else if(names(p_dots$var_exp_spline) == "hladr"){
      var_tbl_add <- dataset_list[[p_dots$dataset_name]]$hladr %>%
        dplyr::mutate(SubjectID = stringr::str_sub(fcs, end = 6),
                      VisitType = purrr::map_chr(fcs, function(fcs_ind){
                        fcs_ind <- stringr::str_sub(fcs_ind, start = 8)
                        first_us_loc <- stringr::str_locate(fcs_ind, "_")[1,"start"][[1]]
                        stringr::str_sub(fcs_ind, end = first_us_loc - 1)
                      })) %>%
        dplyr::select(SubjectID, VisitType, stim, hladr_med_diff)

      names(var_tbl_add) <- c(names(var_tbl_add)[-ncol(var_tbl_add)], names(p_dots$var_exp_spline))

      data_raw %<>%
        inner_join(var_tbl_add,
                   by = c("SubjectID", "VisitType", "stim"))
    }

  }

}

#' @title Scale variables for modelling
scale_var <- function(.data, cols = NULL) {
  if (is.null(cols[1])) {
    return(.data)
  }

  if (length(cols) == 0) return(.data)

  if ("tfmttb" %in% cols && "tfmttb" %in% colnames(.data)) {
    .data <- .data %>%
      dplyr::mutate(origMeantfmttb = mean(tfmttb),
                    tfmttb = tfmttb/1e2)
  }
  if ("DaysSinceEntry" %in% cols && "DaysSinceEntry" %in% colnames(.data)) {
    .data <- .data %>%
      dplyr::mutate(DaysSinceEntry = (DaysSinceEntry - mean(DaysSinceEntry)) /
                      sd(DaysSinceEntry))
  }
  .data
}
