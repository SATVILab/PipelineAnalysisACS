#' @title Pre-process data
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import stringr
#' @import lme4
#' @import glmmTMB
#'
#' @export
preprocess <- function(data_raw, p_dots, dir_proj) {
  data_raw <- data_raw %>%
    tibble::as_tibble()

  # remove NA obs
  data_raw <- data_raw %>%
    dplyr::filter(!is.na(resp))

  # add clinical_data
  # ---------------------
  cols_to_add_vec <- c(
    names(p_dots$var_exp_spline),
    p_dots$var_conf,
    p_dots$var_exp,
    p_dots$var_re,
    p_dots$var_offset
  )

  data_raw <- data_raw %>%
    add_clinical_data_and_filter(
      cols_add = cols_to_add_vec[!grepl("^tc~|^none$", cols_to_add_vec)],
      ttb_min = p_dots$ttb_min,
      ttb_max = p_dots$iter$ttb_max
    )

  # add measurements from other datasets
  # ---------------------
  cols_to_add_vec_tc <- cols_to_add_vec[grepl("^tc~", cols_to_add_vec)] %>%
    stringr::str_remove("^tc~")

  # tuberculomics compendium data
  data_raw <- data_raw %>%
    add_tc_assay_data(
      cols_add = cols_to_add_vec_tc,
      cols_join = c("SubjectID", "VisitType"),
      iter = p_dots$iter
    )

  # transform response
  # ---------------------
  data_raw <- data_raw %>%
    trans(trans = p_dots$trans)

  if (p_dots$family %in% c("bin", "betabin")) {
    data_raw[, "resp"] <- data_raw[, "resp"] / data_raw[, "n_cell"]
  }

  # winsorise, if need be
  if (is.logical(p_dots$iter$wins)) {
    wins <- ifelse(p_dots$iter$wins, "wins_y", "wins")
  } else {
    wins <- p_dots$iter$wins
  }

  data_raw <- winsorise(data_raw = data_raw, wins = wins, p_dots = p_dots)

  # scale measurements, if desired
  # ---------------------
  cols_to_scale <- c(
    "DaysSinceEntry",
    "tfmttb"
  )
  data_raw <- data_raw %>%
    scale_var(cols = cols_to_scale)

  # ===========================


  data_raw
}
