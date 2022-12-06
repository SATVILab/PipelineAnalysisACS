#' @title Pre-process data
#'
#' @import ggplot2
#' @importFrom magrittr |>
#' @import stringr
#' @import lme4
#' @import glmmTMB
#'
#' @export
preprocess <- function(data_raw, iter, p_dots, dir_proj) {
  data_raw <- data_raw |>
    tibble::as_tibble()

  # remove NA obs
  data_raw <- data_raw |>
    dplyr::filter(!is.na(iter$var_conf))

  # add clinical_data
  # ---------------------
  cols_to_add_vec <- c(
    iter[["var_exp_spline"]][[1]]$var |> unlist() |> setNames(NULL),
    iter[["var_conf"]] |> unlist() |> setNames(NULL),
    iter[["var_exp"]] |> unlist() |> setNames(NULL),
    iter[["var_re"]] |> unlist() |> setNames(NULL),
    iter[["var_offset"]] |> unlist() |> setNames(NULL)
  )

  data_raw <- data_raw |>
    add_clinical_data_and_filter(
      cols_add = cols_to_add_vec[!grepl("^tc~|^none$", cols_to_add_vec)],
      ttb_min = iter$ttb_min,
      ttb_max = iter$ttb_max
    )

  # add measurements from other datasets
  # ---------------------
  cols_to_add_vec_tc <- cols_to_add_vec[grepl("^tc~", cols_to_add_vec)] |>
    stringr::str_remove("^tc~")

  # tuberculomics compendium data
  data_raw <- data_raw |>
    add_tc_assay_data(
      cols_add = cols_to_add_vec_tc,
      cols_join = c("SubjectID", "VisitType"),
      iter = iter
    )

  # transform response
  # ---------------------
  data_raw <- data_raw |>
    trans(trans = iter$trans)

  if (iter$family %in% c("bin", "betabin")) {
    data_raw[, "resp"] <- data_raw[, "resp"] / data_raw[, "n_cell"]
  }

  # winsorise, if need be
  if (is.logical(iter$wins)) {
    wins <- ifelse(iter$wins, "wins_y", "wins")
  } else {
    wins <- iter$wins
  }

  data_raw <- winsorise(
    data_raw = data_raw,
    wins = wins,
    p_dots = p_dots,
    iter = iter
  )

  # scale measurements, if desired
  # ---------------------
  cols_to_scale <- c(
    "DaysSinceEntry",
    "tfmttb"
  )
  data_raw <- data_raw |>
    scale_var(cols = cols_to_scale)

  # ===========================


  data_raw
}
