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
  data_raw <- data_raw %>%
    add_clinical_data_and_filter(
      cols_add = c(names(p_dots$var_exp_spline),
                   p_dots$var_conf,
                   p_dots$var_exp,
                   p_dots$var_re,
                   p_dots$var_offset),
      ttb_min = p_dots$ttb_min,
      ttb_max = p_dots$iter$ttb_max
      )

  # choose only progressor or non-prog samples


  # add measurements from other datasets
  # ---------------------
  # TODO: write function to add this

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
