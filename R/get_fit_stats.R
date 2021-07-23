# so, we want (regardless of the model package) the following:
# - the coefficients, their standard errors and their p-values (if output automatically)
# - the wald stats for the explanatory variables (grouped, if need be)
# - the lr stats between models (assumed nested)
# - the kw stats (where you slice continuous variables into a variety of cutpoints,
# perhaps either from the spline (if knots or df provided))
# - mw/kw stats for categorical variables
# - slice into even groups for kw if a continuous variable

#' @export
get_fit_stats <- function(data_raw, data_mod, dir_proj, params_dots, fit_obj){

  # Check if to rerun
  # ==============================

  # check if old output is to be ignored
  ignore_old_output <- switch(as.character(identical(params_dots$ignore_old_output, "")),
                              "TRUE" = FALSE,
                              "FALSE" = any(stringr::str_detect('get_fit|all', params_dots$ignore_old_output)))

  # return if able
  if(!ignore_old_output){
    results_list <- return_rds(file.path(dir_proj, "results"))
    if(length(results_list) > 0) return(results_list)
  }


  # get stats if old ignored or unavailable
  # ==============================

  if(dir.exists(file.path(dir_proj, "results"))){
    unlink(file.path(dir_proj, "results"), recursive = TRUE)
  }

  # prep
  # -------------

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
  var_test_list <- list(params_dots$var_exp, names(params_dots$var_exp_spline),
                       c(params_dots$var_exp, names(params_dots$var_exp_spline)))

  # get Wald stats
  wald_tbl <- modutils::test_wald(fit_obj$full, var_test_list, match_cond = "any")

  # append results
  results_list %<>% append(list('wald' = wald_tbl))

  # lr test
  # -------------
  if(length(setdiff(names(fit_obj), "full")) > 0 && length(fit_obj) > 0){
    var_name_vec <- map_chr(names(fit_obj[-which(names(fit_obj) == 'full')]),
                            function(x){
                              ifelse(x == "null", paste0(c(params_dots$var_exp, names(params_dots$var_exp_spline)),
                                                         collapse = "; "),
                                     x)
                            })
    lr_tbl <- modutils::test_lr(fit_obj$full, fit_obj[-which(names(fit_obj) == 'full')],
                                var_name_vec)

    results_list %<>% append(list('lr' = lr_tbl))
  }

  # save results
  analysispipeline::save_objects(obj_list = results_list,
                                 dir_proj = dir_proj,
                                 dir_sub = "results")

  results_list

}
