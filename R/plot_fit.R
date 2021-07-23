# ====================================

#' @title Fit plots
#'
#' @export
plot_fit <- function(data_raw, data_mod, dir_proj, params_dots, fit_obj,
                        fit_stats){

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
  ignore_old_output <- switch(as.character(identical(params_dots$ignore_old_output, "")),
                              "TRUE" = FALSE,
                              "FALSE" = any(stringr::str_detect('plot_fit|all', params_dots$ignore_old_output)))

  # return if able
  if(!ignore_old_output){
    if(length(list.files(file.path(dir_proj, "results"),
                          pattern = 'png')) > 0) return(invisible(TRUE))
  }

  # =============================
  # Preparation
  # =============================

  # extract model
  mod_full <- fit_obj$full

  var_plot <- c(params_dots$var_conf, names(params_dots$var_exp_spline),
                params_dots$var_exp)
  n_var <- length(var_plot)
  fit_height <- ceiling(n_var / 2) * 8
  fit_width <- ifelse(n_var >= 2, 19, 9.5)

  if(!dir.exists(file.path(dir_proj, "results"))){
    dir.create(file.path(dir_proj, "results"), recursive = TRUE)
  }

  png(filename = file.path(dir_proj, "results", "p_fit_auto.png"), units = 'cm',
      height = fit_height, width = fit_width, res = 300,
      pointsize = 10)

  # add data_mod to global env and ensure removal upon function exit
  # --------------------------
  on.exit(suppressWarnings(try(rm(data_mod, envir = .GlobalEnv), silent = TRUE)))
  if(exists('data_mod', envir = .GlobalEnv)){ # if there is already data_mod
    # in global env, then set data_mod back to that upon exit
    data_mod_old_akx <- get('data_mod', envir = .GlobalEnv)
    on.exit(assign('data_mod', data_mod_old_akx, envir = .GlobalEnv),
            add = TRUE)
  }
  assign(x = "data_mod", value = data_mod, envir = .GlobalEnv)

  # create standard predictorEffects plot
  plot_out <- try(plot(effects::predictorEffects(mod = mod_full)))
  invisible(dev.off())
  if(identical(class(plot_out), 'try-error')){
    rm(file.path(dir_proj, "results", "p_fit_auto.png"))
  }
  rm(plot_out)

  # create plots of var_exp
  # -----------------------------

  if(!is.null(params_dots$var_int) &&
     'Progressor' %in% params_dots$var_int &&
     is.numeric(data_mod[[params_dots$var_int[params_dots$var_int != 'Progressor']]]) &&
     !params_dots$var_int[params_dots$var_int != 'Progressor'] == 'tfmttb'){


    eff_obj <- effects::Effect(params_dots$var_int,
                               mod = mod_full,
                               xlevels = 50)

    eff_tbl <- eff_obj %>%
      as.data.frame() %>%
      as_tibble()

    np_var <- params_dots$var_int[params_dots$var_int != 'Progressor']
    var_alt_ind <- which(colnames(eff_tbl) == np_var)
    colnames(eff_tbl)[var_alt_ind] <- 'var_alt'

    n_cell_ind <- switch(as.character(is.null(params_dots$var_offset)),
                         "TRUE" = FALSE,
                         "FALSE" = stringr::str_detect(params_dots$var_offset, 'n_cell'))

    if(n_cell_ind){
      n_cell <- exp(eff_obj$offset)
      eff_tbl %<>%
        mutate(fit = fit/n_cell*1e2,
               lower = lower/n_cell*1e2,
               upper = upper/n_cell*1e2)
    }

    bounds_tbl <- data_mod %>%
      group_by(Progressor) %>%
      summarise(min = min(!!ensym(np_var)),
                max = max(!!ensym(np_var)))
    min_vec <- setNames(bounds_tbl$min, bounds_tbl$Progressor)
    max_vec <- setNames(bounds_tbl$max, bounds_tbl$Progressor)
    min_prog <- eff_tbl %>%
      filter(Progressor == "yes") %>%
      filter(var_alt <= min_vec["yes"] | var_alt <= min(var_alt)) %>%
      arrange(desc(var_alt)) %>%
      slice(1) %>%
      extract2('var_alt')
    max_prog <- eff_tbl %>%
      filter(Progressor == "yes") %>%
      filter(var_alt >= max_vec["yes"] | var_alt >= max(var_alt)) %>%
      arrange(var_alt) %>%
      slice(1) %>%
      extract2('var_alt')
    min_ctrl <- eff_tbl %>%
      filter(Progressor == "no") %>%
      filter(var_alt <= min_vec["no"] | var_alt <= min(var_alt)) %>%
      arrange(desc(var_alt)) %>%
      slice(1) %>%
      extract2('var_alt')
    max_ctrl <- eff_tbl %>%
      filter(Progressor == "no") %>%
      filter(var_alt >= max_vec["no"]  | var_alt >= max(var_alt)) %>%
      arrange(var_alt) %>%
      slice(1) %>%
      extract2('var_alt')
    min_vec <- c("yes" = min_prog, 'no' = min_ctrl)
    max_vec <- c("yes" = max_prog, 'no' = max_ctrl)

    plot_tbl_prog <- eff_tbl %>%
      filter(Progressor == 'yes') %>%
      filter(var_alt >= min_vec['yes'],
             var_alt <= max_vec['yes'])
    plot_tbl_ctrl <- eff_tbl %>%
      filter(Progressor == 'no') %>%
      filter(var_alt >= min_vec['no'],
             var_alt <= max_vec['no'])
    plot_tbl <- plot_tbl_prog %>%
      bind_rows(plot_tbl_ctrl)

    p <- ggplot(data = plot_tbl,
                mapping = aes(x = var_alt, col = Progressor, fill = Progressor)) +
      background_grid(major = 'xy', minor = 'y') +
      #scale_x_reverse() +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.5) +
      geom_line(aes(y = fit)) +
      scale_fill_manual(values = c("yes" = "orange",
                                   "no" = "dodgerblue")) +
      scale_colour_manual(values = c("yes" = "orange",
                                     "no" = "dodgerblue")) +
      labs(x = np_var, y = "Response")

    lr_tbl <- readRDS(file.path(dir_proj, "results", "lr.rds"))
    lr_tbl %<>%
      mutate(var = map_chr(var, function(x){
        x %>%
          stringr::str_replace("; ", " and/or ")
      })) #%>%
      #arrange(var)

    res_text <- paste0(lr_tbl$var, ": ", signif(lr_tbl$p, digits = 3))

    data_raw_plot <- data_mod %>%
      mutate(ttb =  params_dots$max_ttb - (tfmttb * 1e2)) %>%
      mutate(resp = purrr::map_dbl(1:nrow(data_mod), function(i){
        if(!n_cell_ind) return(data_mod$resp[i])
        data_mod$resp[i]/data_mod$n_cell[i]*1e2
      })) %>%
      select(!!np_var, resp, Progressor, ttb) %>%
      mutate(ttb_cat = ifelse(ttb > 360, "360-", ""),
             ttb_cat = ifelse(ttb > 180 & ttb <= 360, "180-360", ttb_cat),
             ttb_cat = ifelse(ttb <= 180, "0-180", ttb_cat)) %>%
      rename(`Time to TB` = ttb_cat)
    var_alt_ind <- which(colnames(data_raw_plot) == np_var)
    colnames(data_raw_plot)[var_alt_ind] <- 'var_alt'

    p_raw <- p +
      geom_point(data = data_raw_plot,
                 aes(y = resp, shape = `Time to TB`, alpha = Progressor,
                     size = Progressor)) + #,
                 #size = 1.25) +
      scale_shape_manual(values = c("360-" = 6,
                                    "180-360" = 5,
                                    "0-180" = 12)) +
      scale_alpha_manual(values = c('no' = 0.33, 'yes' = 1)) +
      scale_size_manual(values = c('no' = 1, 'yes' = 1.6))

    if(n_cell_ind){
      p <- p + coord_cartesian(ylim = c(0, min(100, max(eff_tbl$upper))))
      p_raw <- p_raw + coord_cartesian(ylim = c(0,
                                                min(100, max(eff_tbl$upper,
                                                             data_raw_plot$resp))))
    }

    p_list_res <- map(list(p, p_raw), function(p){
      p <- ggdraw(p)
      for(i in seq_along(res_text)){
        p <- p +
          draw_text(res_text[i], x = 0.1675, y = 0.95 - 0.05 * (i-1), size = 12 , hjust = 0)
      }
      p
    }) %>%
      setNames(c("p_fit_manual_without_raw_data",
                 "p_fit_manual_with_raw_data"))


    analysispipeline::save_objects(obj_list = p_list_res[2],
                                   dir_proj = dir_proj,
                                   dir_sub = 'results',
                                   empty = FALSE,
                                   width = 19,
                                   height = 15)


  }



  if(identical(params_dots$var_exp, "Progressor") && identical(names(params_dots$var_exp_spline), 'tfmttb')){


    eff_obj <- effects::Effect(c("tfmttb", "Progressor"),
                               mod = mod_full,
                               xlevels = 50)

    eff_tbl <- eff_obj %>%
      as.data.frame() %>%
      as_tibble()

    eff_tbl %<>%
      filter(Progressor == 'yes' | (Progressor == 'no' & tfmttb == 0))

    n_cell_ind <- switch(as.character(is.null(params_dots$var_offset)),
                         "TRUE" = FALSE,
                         "FALSE" = stringr::str_detect(params_dots$var_offset, 'n_cell'))

    if(n_cell_ind){
      n_cell <- exp(eff_obj$offset)
      eff_tbl %<>%
        mutate(fit = fit/n_cell*1e2,
               lower = lower/n_cell*1e2,
               upper = upper/n_cell*1e2)
    }

    eff_tbl_prog <- eff_tbl %>%
      filter(Progressor == 'yes')

    eff_tbl_ctrl <- map_df(eff_tbl_prog$tfmttb, function(x){
      eff_tbl %>%
        filter(Progressor == 'no',
               tfmttb == 0) %>%
        mutate(tfmttb = .env$x)
    })

    eff_tbl_plot <- eff_tbl_prog %>%
      bind_rows(eff_tbl_ctrl)

    p <- ggplot(data = eff_tbl_plot %>%
                  mutate(ttb = params_dots$max_ttb - (tfmttb * 1e2)),
                mapping = aes(x = ttb, col = Progressor, fill = Progressor)) +
      background_grid(major = 'xy', minor = 'y') +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.5) +
      geom_line(aes(y = fit)) +
      scale_x_reverse() +
      scale_fill_manual(values = c("yes" = "orange",
                                   "no" = "dodgerblue")) +
      scale_colour_manual(values = c("yes" = "orange",
                                     "no" = "dodgerblue")) +
      labs(x = "Days to TB", y = "Response")

    lr_tbl <- readRDS(file.path(dir_proj, "results", "lr.rds"))
    var_lab_vec <- c("Progressor; tfmttb" = "Progressor and/or time",
                     "Progressor" = "Progressor",
                     "tfmttb" = "Time")
    lr_tbl %<>%
      mutate(var = var_lab_vec[var],
             var = factor(.data$var, levels = var_lab_vec)) %>%
      arrange(var)

    res_text <- paste0(lr_tbl$var, ": ", signif(lr_tbl$p, digits = 3))


    data_raw_plot <- data_mod %>%
      mutate(resp = purrr::map_dbl(1:nrow(data_mod), function(i){
        if(!n_cell_ind) return(data_mod$resp[i])
        data_mod$resp[i]/data_mod$n_cell[i]*1e2
        }),
        ttb = params_dots$max_ttb - (tfmttb * 1e2)) %>%
      select(ttb, resp, Progressor)

    p_raw <- p +
      geom_point(data = data_raw_plot,
                 aes(y = resp))

    if(n_cell_ind){
      p <- p + coord_cartesian(ylim = c(0, min(100, max(eff_tbl_plot$upper))))
      p_raw <- p_raw + coord_cartesian(ylim = c(0, min(100, max(eff_tbl_plot$upper, data_raw_plot$resp))))
    }

    p_list_res <- map(list(p, p_raw), function(p){
      ggdraw(p) +
        draw_text(res_text[1], x = 0.1675, y = 0.95, size = 12 , hjust = 0) +
        draw_text(res_text[2], x = 0.1675, y = 0.95 - 0.05, size = 12, hjust = 0) +
        draw_text(res_text[3], x = 0.1675, y = 0.95 - 2 * 0.05, size = 12, hjust = 0)
    }) %>%
      setNames(c("p_fit_manual_without_raw_data",
                 "p_fit_manual_with_raw_data"))


    analysispipeline::save_objects(obj_list = p_list_res,
                                   dir_proj = dir_proj,
                                   dir_sub = 'results',
                                   empty = FALSE,
                                   width = 19,
                                   height = 15)
    analysispipeline::save_objects(obj_list = p_list_res,
                                   dir_proj = dir_proj,
                                   dir_sub = 'results',
                                   empty = FALSE,
                                   width = 19,
                                   gg_device = "pdf",
                                   height = 15)


   }

  return(invisible(TRUE))

}
