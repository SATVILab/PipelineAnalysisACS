# ====================================

#' @title Fit plots
#'
#' @export
#'
#' @importFrom rlang !! ensym :=
#' @importFrom splines ns
display <- function(data_raw, data_mod, dir_proj,
                    iter,
                    p_dots, fit_obj,
                    fit_stats) {
  if (identical(class(fit_obj$full), "try-error")) {
    return(invisible(TRUE))
  }
  fit_obj <- fit_obj[
    purrr::map_lgl(
      fit_obj,
      function(x) !identical(class(x), "try-error")
    )
  ]

  ggplot2::theme_set(cowplot::theme_cowplot())

  library(splines)
  on.exit(suppressWarnings(try(
    detach("package:splines", unload = TRUE),
    silent = TRUE
  )))


  # goal is to plot residuals and CI's for every "point"

  # save auto effects pkg plot
  # -----------------
  plot_display_auto(
    mod = fit_obj$full,
    data_mod = data_mod,
    var = c(
      iter$var_conf, names(iter$var_exp_spline),
      iter$var_exp
    ),
    dir_save = iter$dir_stg
  )

  # create plots of var_exp
  # -----------------------------

  if ("cd4_ifng_freq" %in% names(iter$var_exp_spline)) {
    eff_obj <- effects::Effect(c("cd4_ifng_freq"),
      mod = mod_full,
      xlevels = 50
    )
  }

  # plot response against interaction with progressor
  # --------------

  is_int <- !is.null(iter$var_int[[1]])
  var_int_non_p <- setdiff(iter$var_int[[1]], "Progressor")
  if (is_int) {
    is_int_prog <- "Progressor" %in% iter$var_int[[1]]
    # prog included
    is_int_prog_n_num <- is.numeric(data_mod[[var_int_non_p]]) # numeric
    plot_int <- is_int && is_int_prog && is_int_prog_n_num
  } else {
    plot_int <- FALSE
    var_int_non_p <- "tfmttb"
  }

  is_tfmttb_mod <- grepl(
    "tfmttb",
    names(iter$var_exp_spline),
  )
  if (plot_int || is_tfmttb_mod) {
    axis_lab_x <- ifelse(var_int_non_p == "tfmttb",
      "Days to TB",
      stringr::str_to_upper(var_int_non_p)
    )
    axis_lab <- c(axis_lab_x, iter$var_dep)

    purrr::walk(c(Inf, 5), function(max_sd) {
      iter <- rlang::caller_env(3)$iter
      p_list <- plot_disp_int_cat_num(
        mod = fit_obj$full,
        .data = data_mod,
        data_nm = "data_mod",
        var_num = var_int_non_p,
        var_cat = "Progressor",
        var_offset = switch(
          as.character(iter[["var_offset"]] == "none"),
          "TRUE" = NULL,
          iter[["var_offset"]],
        var_dep = iter[["var_dep"]],
        max_sd = max_sd,
        cat_to_col = c(
          "yes" = "orange",
          "no" = "dodgerblue"
        ),
        axis_lab = axis_lab,
        axis_x_reverse = var_int_non_p == "tfmttb",
        add_test = "lr",
        dir_test = file.path(
          dirname(iter$dir_stg),
          "extr"
        )
      )

      p_list <- p_list[grepl("raw", names(p_list))]

      names(p_list) <- paste0(
        names(p_list),
        ifelse(max_sd == Inf, "", paste0("_maxsd", max_sd))
      )

      purrr::walk(c("pdf", "png"), function(gd) {
        pipeline::save_objects(
          obj_list = p_list,
          dir_proj = iter$dir_stg,
          empty = FALSE,
          width = 19,
          height = 15,
          gg_device = gd
        )
      })
    })
  }
  return(invisible(TRUE))
}
