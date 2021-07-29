#' @export
validate <- function(data_raw, data_mod, dir_proj, p_dots, fit_obj){

  theme_set(cowplot::theme_cowplot())

  # ===========================
  # Plot validation plots
  # ===========================

  p_dots <- remove_tc_assay_from_exp_s(p_dots)

  # directory to save to
  dir_save <- p_dots$dir_stg

  # object to save plots to
  p_list <- list()

  # gather variables to make residuals against
  var_plot_exp_vec <- c(p_dots$var_exp, names(p_dots$var_exp_spline))
  var_plot_vec <- c(".no_exp_var", var_plot_exp_vec) # adds variable meaning "no variables"

  # get residuals
  var_extract_vec <- c(var_plot_exp_vec, p_dots$var_dep)
  data_plot <- data_mod[,var_extract_vec] %>%
    dplyr::bind_cols(modutils::get_resid(fit_obj$full))

  # get fitted values
  data_plot <- data_plot %>%
    dplyr::bind_cols(modutils::get_pred(fit_obj$full)[,".pred_resp"])

  # make univariate residual plots
  # ---------------------------
  p_resid_list <- purrr::map(var_plot_vec, function(var){

    if(var == ".no_exp_var"){
      # plot against N(0,1) dbn for clarity
      dens_res_std <- density(data_plot[[".resid_std"]])
      data_res_std_dens <- tibble::tibble(x = dens_res_std$x, y = dens_res_std$y,
                                  type = "Std residuals")
      data_norm <- tibble::tibble(x = seq(min(data_res_std_dens$x),
                                  max(data_res_std_dens$x),
                                  length.out = nrow(data_res_std_dens)),
                          type = "N(0,1)")
      data_norm$y <- dnorm(data_norm$x)
      data_plot_no_exp <- data_res_std_dens %>%
        dplyr::bind_rows(data_norm)
      return(list(p = ggplot(data_plot_no_exp) +
               cowplot::background_grid(major = 'x') +
               geom_line(aes(x = x, y = y, col = type),
                         show.legend = FALSE) +
               labs(x = "Standardised residuals", y = "Density"),
               "save" = c(height = 10, width = 18)))
    }
    p <- ggplot(data_plot, aes(x = !!ensym(var), y = .resid_std)) +
      geom_hline(yintercept = 0)
    p <- switch(as.character(is.character(data_plot[[var]]) || is.factor(data_plot[[var]])),
                "TRUE" = p + geom_boxplot() +
                  cowplot::background_grid(major = 'y') +
                  ggforce::geom_sina(),
                "FALSE" = p +
                  cowplot::background_grid(major = 'xy') +
                  geom_point() +
                  geom_smooth(se = FALSE)) +
      labs(y = "Standardised residuals")

    list(p = p, "save" = c("height" = 10, width = 18))
  }) %>%
    setNames(paste0("resid-", var_plot_vec))

  purrr::walk(seq_along(p_resid_list), function(i){
    p_save <- setNames(list(p_resid_list[[i]]$p), names(p_resid_list)[i])
    pipeline::save_objects(obj_list = p_save,
                                   dir_proj = dir_save,
                                   height = p_resid_list[[i]]$save["height"],
                                   width = p_resid_list[[i]]$save["width"],
                                   empty = FALSE)
  })


  # make plots of residuals against fitted values
  p_resid_vs_fitted <- purrr::map(var_plot_vec, function(var){
    if(is.numeric(data_plot[[var]])){
      data_plot[[var]] <- cut(data_plot[[var]], breaks = max(min(5, nrow(data_plot[[var]] - 2)), 1))
    }
    p <- ggplot(data_plot, aes(x = .pred_resp,  y = .resid_std)) +
      cowplot::background_grid(major = 'xy') +
      geom_point() +
      geom_smooth() +
      labs(x = "Fitted value", y = "Standardised residuals")
    if(var != '.no_exp_var'){
      if(is.character(data_plot[[var]])){
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ", ncol = 2)")))
        p <- p + labs(title = var)
      } else{
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ", ncol = 2)")))
        p <- p + labs(title = var)
      }
    }
    # save plots with height dependent on number of variables included
    list(p = p, save = c("height" = ifelse(var == '.no_exp_var',
                                           10, 10 * length(unique(data_plot[[var]]))/2), "width" = 18))
  }) %>%
    setNames(paste0("resid_vs_fitted-", var_plot_vec))

  purrr::walk(seq_along(p_resid_vs_fitted), function(i){
    p_save <- setNames(list(p_resid_vs_fitted[[i]]$p), names(p_resid_vs_fitted)[i])
    pipeline::save_objects(obj_list = p_save,
                                   dir_proj = dir_save,
                                   height = p_resid_vs_fitted[[i]]$save["height"],
                                   width = p_resid_vs_fitted[[i]]$save["width"],
                                   empty = FALSE)
  })

  # make plots of residuals against actual
  p_resid_vs_actual <- purrr::map(var_plot_vec, function(var){

    if(is.numeric(data_plot[[var]])){
      data_plot[[var]] <- cut(data_plot[[var]], breaks = max(min(5, nrow(data_plot[[var]] - 2)), 1))
    }
    p <- ggplot(data_plot, aes(x = !!sym(p_dots$var_dep),  y = .resid_std)) +
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'xy') +
      geom_point() +
      geom_smooth(se = FALSE) +
      labs(x = "Fitted value", y = "Standardised residuals")
    if(var != '.no_exp_var'){
      if(is.character(data_plot[[var]])){
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ")")))
        p <- p + labs(title = var)
      } else{
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ")")))
        p <- p + labs(title = var)
      }
    }
    # save plots with height dependent on number of variables included
    list(p = p, save = c("height" = ifelse(var == '.no_exp_var',
                                           10, 10 * length(unique(data_plot[[var]]))/2), "width" = 18))
  }) %>%
    setNames(paste0("resid_vs_actual-", var_plot_vec))

  purrr::walk(seq_along(p_resid_vs_actual), function(i){
    p_save <- setNames(list(p_resid_vs_actual[[i]]$p), names(p_resid_vs_actual)[i])
    pipeline::save_objects(obj_list = p_save,
                                   dir_proj = dir_save,
                                   height = p_resid_vs_actual[[i]]$save["height"],
                                   width = p_resid_vs_actual[[i]]$save["width"],
                                   empty = FALSE)
  })

  # make plots of fitted vs actual values
  p_fitted_vs_actual <- purrr::map(var_plot_vec, function(var){
    if(is.numeric(data_plot[[var]])){
      data_plot[[var]] <- cut(data_plot[[var]], breaks = max(min(5, nrow(data_plot[[var]] - 2)), 1))
    }
    p <- ggplot(data_plot, aes(x =  !!sym(p_dots$var_dep),  y = .pred_resp)) +
      geom_abline(intercept = 0, slope = 1) +
      cowplot::background_grid(major = 'xy') +
      geom_point() +
      geom_smooth(se = FALSE) +
      labs(x = "Actual value", y = "Fitted value")
    if(var != '.no_exp_var'){
      if(is.character(data_plot[[var]])){
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ")")))
        p <- p + labs(title = var)
      } else{
        p <- eval(parse(text = paste0("p + facet_wrap(~", var, ")")))
        p <- p + labs(title = var)
      }
    }
    # save plots with height dependent on number of variables included
    list(p = p, save = c("height" = ifelse(var == '.no_exp_var',
                                           10, 10 * length(unique(data_plot[[var]]))/2), "width" = 18))
  }) %>%
    setNames(paste0("fitted_vs_actual-", var_plot_vec))

  purrr::walk(seq_along(p_fitted_vs_actual), function(i){
    p_save <- setNames(list(p_fitted_vs_actual[[i]]$p), names(p_fitted_vs_actual)[i])
    pipeline::save_objects(obj_list = p_save,
                                   dir_proj = dir_save,
                                   height = p_fitted_vs_actual[[i]]$save["height"],
                                   width = p_fitted_vs_actual[[i]]$save["width"],
                                   empty = FALSE)
  })

  # =======================
  # get non-parametric tests
  # =======================

  var_exp_vec <- setdiff(var_plot_vec, ".no_exp_var")
  var_exp_list <- as.list(var_exp_vec) %>% append(list(var_exp_vec))
  var_spline_vec <- names(p_dots$var_exp_spline)
  knot_list <- switch(as.character(is.null(p_dots$var_exp_spline)),
                      "TRUE" = list(),
                      "FALSE" = map(p_dots$var_exp_spline, function(x){
                        if('knots' %in% names(x$params)) return(x$params$knots)
                        NULL
                      }) %>%
                        setNames(names(p_dots$var_exp_spline)))

  break_list <- map(var_exp_vec, function(var){
    if(is.character(data_mod[[var]]) || is.logical(data_mod[[var]])) return(unique(data_mod[[var]]))
    if(is.factor(data_mod[[var]])) return(levels(data_mod[[var]]))
    if(var %in% var_spline_vec){
      if('knots' %in% names(p_dots$var_exp_spline[[var]]$params)){
        break_vec <- c(min(data_mod[[var]]), p_dots$var_exp_spline[[var]]$params$knots,
                       max(data_mod[[var]]))
        break_vec <- unique(break_vec)
        range_len <- diff(range(break_vec))
        break_vec[1] <- break_vec[1] - 0.01 * range_len
        break_vec[length(break_vec)] <- break_vec[length(break_vec)] + 0.01 * range_len
        return(break_vec)
      }
    }
    quant_vec <- seq(0,1, length.out = min(5, length(unique(data_mod[[var]]))))
    break_vec <- quantile(data_mod[[var]], quant_vec)
    range_len <- diff(range(break_vec))
    break_vec[1] <- break_vec[1] - 0.01 * range_len
    break_vec[length(break_vec)] <- break_vec[length(break_vec)] + 0.01 * range_len
    break_vec

    }) %>%
    setNames(var_exp_vec)
  #browser()
  kw_obj <- map_df(var_exp_list, function(var_exp){
    #print(var_exp)
    resp_vec <- data_mod[[p_dots$var_dep]]
    if(!is.null(p_dots$var_offset)){
      resp_vec <- resp_vec/data_mod[[p_dots$var_offset]]
    }
    for(i in seq_along(var_exp)){
      if(i == 1){
        grp_vec <- data_mod[[var_exp[i]]]
        if(is.logical(grp_vec)) grp_vec <- as.character(grp_vec)
        if(is.numeric(grp_vec)){
          grp_vec <- cut(grp_vec, breaks = break_list[[var_exp[i]]])
        }
      } else{
        grp_vec_curr <- data_mod[[var_exp[i]]]
        if(is.logical(grp_vec_curr)) grp_vec_curr <- as.character(grp_vec_curr)
        if(is.numeric(grp_vec_curr)){
          grp_vec_curr <- cut(grp_vec_curr, breaks = break_list[[var_exp[i]]])
        }
        grp_vec <- paste0(grp_vec, "_", grp_vec_curr)
      }
    }
    try(kw_obj <- kruskal.test(resp_vec, grp_vec), silent = TRUE)
    var_print <- paste0(var_exp, collapse = "; ")

    out_tbl <- try(tibble::tibble(var = var_print, stat = kw_obj$statistic, df = length(unique(grp_vec)) - 1,
           p = kw_obj$p.value))
    if(class(out_tbl) == 'try-error') return(tibble(var = NA_character_, stat = NA_real_, df = NA_integer_))
    out_tbl
  })

  pipeline::save_objects(obj_list = list(kw = kw_obj),
                                 dir_proj = dir_save,
                                 empty = FALSE)
  invisible(TRUE)

}
