# so, we want (regardless of the model package) the following:
# - the coefficients, their standard errors and their p-values (if output automatically)
# - the wald stats for the explanatory variables (grouped, if need be)
# - the lr stats between models (assumed nested)
# - the kw stats (where you slice continuous variables into a variety of cutpoints,
# perhaps either from the spline (if knots or df provided))
# - mw/kw stats for categorical variables
# - slice into even groups for kw if a continuous variable

#' @export
extract <- function(data_raw, data_mod, dir_proj, p_dots, fit_obj){

  if (identical(class(fit_obj$full), "try-error")) {
    return(invisible(TRUE))
  }
  fit_obj <- fit_obj[
    purrr::map_lgl(
      fit_obj,
      function(x) !identical(class(x), "try-error"))
  ]

  # prep
  # -------------

  p_dots <- remove_tc_assay_from_exp_s(p_dots)

  # list to save results to
  results_list <- list()

  # coefficients
  # -------------
  coef_tbl <- modutils::get_coef(fit_obj$full)
  # append results
  results_list %<>% append(list('coef' = coef_tbl))

  # wald
  # -------------

  # variables to test for in wald test
  var_test_list <- list(p_dots$var_exp, names(p_dots$var_exp_spline),
                        c(p_dots$var_exp, names(p_dots$var_exp_spline)))
  var_test_list <- var_test_list[!vapply(var_test_list, is.null, logical(1))]

  # get Wald stats
  wald_tbl <- modutils::test_wald(
    fit_obj$full, var_test_list, match_cond = "any"
  )

  # append results
  results_list %<>% append(list('wald' = wald_tbl))

  # lr test
  # -------------
  if(length(setdiff(names(fit_obj), "full")) > 0 && length(fit_obj) > 0){
    var_name_vec <- purrr::map_chr(
      names(fit_obj[-which(names(fit_obj) == 'full')]),
      function(x){
        ifelse(
          x == "null",
          paste0(c(p_dots$var_exp, names(p_dots$var_exp_spline)),
                 collapse = "; "),
          x)
      })
    lr_tbl <- modutils::test_lr(
      fit_obj$full,
      fit_obj[-which(names(fit_obj) == 'full')],
      var_name_vec)

    results_list %<>% append(list('lr' = lr_tbl))
  }

  # save results
  pipeline::save_objects(obj_list = results_list,
                         dir_proj = p_dots$dir_stg,
                         empty = FALSE)
  pipeline::save_objects(obj_list = list("fit_stats" = results_list),
                         dir_proj = p_dots$dir_stg,
                         empty = FALSE)

  results_list

}
