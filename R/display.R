# ====================================

#' @title Fit plots
#'
#' @export
#'
#' @importFrom rlang !! ensym :=
#' @importFrom splines ns
display <- function(data_raw, data_mod, dir_proj,
                    p_dots, fit_obj,
                    fit_stats){

  if (identical(class(fit_obj$full), "try-error")) {
    return(invisible(TRUE))
  }
  fit_obj <- fit_obj[
    purrr::map_lgl(
      fit_obj,
      function(x) !identical(class(x), "try-error"))
  ]

  p_dots <- remove_tc_assay_from_exp_s(p_dots)

  theme_set(cowplot::theme_cowplot())

  library(splines)
  on.exit(suppressWarnings(try(
    detach("package:splines", unload = TRUE),
    silent = TRUE
  )))


  #return(invisible(TRUE))

  # goal is to plot residuals and CI's for every
  # not sure what "every" is.
  # I suppose I mean that we should plot them for every level
  # of every factor/character var in c(var_exp, names(var_exp_spline)) and
  # across the range of all continuous variables,
  # and for the combination of the two:
  # for both factors, then for each crossed combination.
  # for one factor and the other continuous, then the response plotted for the
  # continuous var .
  # - could possibly log-transform automatically if the data are VERY skew.
  # or the bulk of the data are very small relative to the total.
  #
  # should always include the intercept, definitely.

  # check if old output is to be ignored

  # save auto effects pkg plot
  # -----------------
  plot_display_auto(
    mod = fit_obj$full,
    data_mod = data_mod,
    var = c(p_dots$var_conf, names(p_dots$var_exp_spline),
            p_dots$var_exp),
    dir_save = p_dots$dir_stg
  )




  # create plots of var_exp
  # -----------------------------

  if ("cd4_ifng_freq" %in% names(p_dots$iter$var_exp_spline)) {
    eff_obj <- effects::Effect(c("cd4_ifng_freq"),
                               mod = mod_full,
                               xlevels = 50)
  }

  # plot response against interaction with progressor
  # --------------

  is_int <- !is.null(p_dots$var_int) # interaction present
  var_int_non_p <- setdiff(p_dots$var_int, "Progressor")
  if (is_int) {
    is_int_prog <- "Progressor" %in% p_dots$var_int
    # prog included
    is_int_prog_n_num <- is.numeric(data_mod[[var_int_non_p]]) # numeric
    plot_int <- is_int && is_int_prog && is_int_prog_n_num
  } else {
    plot_int <- FALSE
    var_int_non_p <- "tfmttb"
  }

  is_tfmttb_mod <- identical(names(p_dots$var_exp_spline),
                             "tfmttb")
  if(plot_int || is_tfmttb_mod){
    axis_lab_x <- ifelse(var_int_non_p == "tfmttb",
                         "Days to TB",
                         stringr::str_to_upper(var_int_non_p))
    axis_lab <- c(axis_lab_x, p_dots$var_dep)
    p_list <- plot_disp_int_cat_num(
      mod = fit_obj$full,
      .data = data_mod,
      data_nm = "data_mod",
      var_num = var_int_non_p,
      var_cat = "Progressor",
      var_offset = p_dots$var_offset,
      var_dep = p_dots$var_dep,
      cat_to_col = c("yes" = "orange",
                     "no" = "dodgerblue"),
      axis_lab = axis_lab,
      axis_x_reverse = var_int_non_p == "tfmttb",
      add_test = "lr",
      dir_test = file.path(dirname(p_dots$dir_stg),
                           "extr"))

    purrr::walk(c("pdf", "png"), function(gd) {
      pipeline::save_objects(obj_list = p_list,
                             dir_proj = p_dots$dir_stg,
                             empty = FALSE,
                             width = 19,
                             height = 15,
                             gg_device = gd)
    })
  }
  return(invisible(TRUE))
}
