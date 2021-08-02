#' @export
filter_using_list <- function(.data, filter_list) {
  for (i in seq_along(filter_list)) {
    nm <- names(filter_list)[i]
    value <- filter_list[[i]]
    .data <- .data[.data[[nm]] %in% value, ]
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
                    count_stim, count_uns) %>%
      dplyr::filter(grepl("\\+", cyt_combn)),
    stop(paste0(
        paste0(cyt_response_type, collapse = "/"),
        " not available for cyt_response_type in prep_bs_freq")
        )
    )



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
  ttb_max <- ifelse(is.na(ttb_max),
                    max(clinical_data_add$timeToTBFromVisit,
                        na.rm = TRUE),
                    ttb_max)
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
                              cols_add,
                              cols_join) {

  cn_orig <- colnames(.data)

  if (length(cols_add) == 0) return(.data)

  cols_add_split <- stringr::str_split(cols_add, "~")
  assay_vec <- purrr::map_chr(cols_add_split, function(x) x[1]) %>%
    unique
  assay_to_cols_add <- purrr::map(assay_vec, function(x) {
    purrr::map_chr(cols_add_split, function(x) x[2])[
      purrr::map_chr(cols_add_split, function(x) x[1]) == x
    ]
  }) %>%
    setNames(assay_vec)

  for (assay in assay_vec) {
    .add_tc_assay_data <- switch(
      assay,
      "soma" = .add_tc_assay_data_soma,
      "risk6" = .add_tc_assay_data_risk6,
      "flow_ifng" = .add_tc_assay_data_flow_ifng,
      stop("assay not recognised")
    )
    .data <- .add_tc_assay_data(
      .data = .data,
      cols_add = assay_to_cols_add[[assay]],
      cols_join = cols_join
      )
  }

  .data[, c(
    cols_join,
    setdiff(colnames(.data), c(cn_orig)),
    setdiff(cn_orig, cols_join))]

}

.add_tc_assay_data_flow_ifng <- function(.data,
                                    cols_add,
                                    cols_join) {
  var_tbl_add <-
    TuberculomicsCompendium::tidy_wb_cd4p_infgp_cells_hladr_mfi

  var_tbl_add[,cols_add] <- var_tbl_add$CD4IFNg / var_tbl_add$CD4 * 1e2
  var_tbl_add <- var_tbl_add[,unique(c(cols_join, cols_add))]
  .data %>%
    dplyr::inner_join(
      var_tbl_add,
      by = cols_join)
}

.add_tc_assay_data_soma <- function(.data,
                                    cols_add,
                                    cols_join) {

  var_tbl_add <- TuberculomicsCompendium::soma_data_tidy %>%
    dplyr::filter(Soma_Target %in% cols_add) %>%
    dplyr::group_by_at(c(cols_join, "Soma_Target")) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  var_tbl_add <- var_tbl_add %>%
    tidyr::pivot_wider(
      id_cols = cols_join,
      names_from = "Soma_Target",
      values_from = "Soma_TransformedReadout"
    )

  var_tbl_add <- var_tbl_add[, c(cols_join, cols_add)]
  new_col_vec <- setdiff(colnames(var_tbl_add), cols_join)

  .data <- .data %>%
    dplyr::inner_join(
      var_tbl_add,
      by = cols_join)

  .data[, c(
    cols_join,
    new_col_vec,
    setdiff(colnames(.data),
            c(cols_join, new_col_vec)))]
}

.add_tc_assay_data_risk6 <- function(.data,
                                     cols_add,
                                     cols_join) {
  .data <- .data %>%
    dplyr::inner_join(
      TuberculomicsCompendium::signature_6gene,
      by = cols_join
    ) %>%
    dplyr::rename(risk6 = sig6gene_CorScore)
  .data[,c(cols_join, setdiff(colnames(.data), cols_join))]
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


remove_tc_assay_from_exp_s <- function(p_dots) {
  nm_vec <- names(p_dots$var_exp_spline)
  for (i in seq_along(nm_vec)) {
    nm <- nm_vec[i]
    if(!grepl("^tc~", nm)) next
    nm_vec[i] <- stringr::str_split(nm, "~")[[1]][3]
  }
  names(p_dots$var_exp_spline) <- nm_vec
  p_dots
}
