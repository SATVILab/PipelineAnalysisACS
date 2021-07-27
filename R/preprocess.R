#' @title Pre-process data
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import stringr
#' @import lme4
#' @import glmmTMB
#'
#' @export
preprocess <- function(data_raw, p_dots, dir_proj){


  data_raw %<>% tibble::as_tibble()

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

  # choose only progressor or non-prog samples


  # add measurements from other datasets
  # ---------------------
  cols_to_add_vec_tc <- cols_to_add_vec[grepl("^tc~", cols_to_add_vec)] %>%
    stringr::str_remove("^tc~")
  data_raw <- data_raw %>%
    add_tc_assay_data(
      cols_add = cols_to_add_vec_tc,
      cols_join = c("SubjectID", "VisitType")
    )
  nm_vec <- names(p_dots$var_exp_spline)
  for (i in seq_along(nm_vec)) {
    nm <- nm_vec[i]
    if(!grepl("^tc~", nm)) next
    nm_vec[i] <- stringr::str_split(nm, "~")[[1]][3]
  }
  names(p_dots$var_exp_spline) <- nm_vec
  assign("p_dots", p_dots, envir = rlang::caller_env())


  # transform response
  # ---------------------
  data_raw <- data_raw %>%
    trans(trans = p_dots$trans)

  # scale measurements, if desired
  # ---------------------
  cols_to_scale <- c("DaysSinceEntry",
                     "tfmttb")
  data_raw <- data_raw %>%
    scale_var(cols = cols_to_scale)

  # ===========================


  data_raw
}
