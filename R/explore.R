# ====================================

#' @title Make exploratory plots
#' @export
explore <- function(data_raw, data_mod, dir_proj, p_dots){


  theme_set(cowplot::theme_cowplot())

  # =============================
  # Preparation
  # =============================


  # miscellaneous
  # -------------------
  # get max_ttb
  max_ttb <- get_max_ttb(data_mod = data_mod, data_raw = data_raw,
                         p_dots = p_dots)

  # get indicator if proportions should be plotted
  prop_ind <- p_dots$response_type == 'count' && p_dots$assay %in% c("cytof", "cytof_ag_flowsom")

  # font size
  y_axis_title_font_size <- 4
  cowplot_theme_font_size <- 12

  # y-axis labels
  if(p_dots$assay %in% c('cytof', 'cytof-faust', 'cytof_flowsom')){
    if(p_dots$response_type == 'count'){
      y_lab <- ifelse(str_detect(p_dots$assay, 'faust'), 'Frequency',
                      paste0(p_dots$stim, "-induced cytokine+ frequency"))
    } else if(p_dots$response_type == 'compass'){
      y_lab <- paste0("Probability of response to ", p_dots$stim)
    } else if(p_dots$response_type == 'hladr_med_diff'){
      y_lab <- paste0("Median difference in HLADR expression level of", p_dots$stim, "+ cells")
    }
  } else y_lab <- paste0(p_dots$resp_type, " response")

  # data for plot
  # -------------------

  data_plot <- data_mod# %>%
    #left_join(data_raw %>%
     #           select(pop, fcs, resp_type, Progressor, SampleID,
     #                  timeToTBFromVisit, timeToTB),
       #       by = c("pop", "fcs", "resp_type", "Progressor", "SampleID"))
  data_plot %<>% filter(!is.na(resp))
  if(prop_ind) data_plot %<>% mutate(resp = resp/n_cell * 100)

  # ==============================
  # Plot - FlowSOM marker plots
  # ==============================



  if(p_dots$assay == 'cytof_ag_flowsom'){

    if(!dir.exists(file.path(dir_proj, 'exp'))) dir.create(file.path(dir_proj, 'exp'), recursive = TRUE)

    path_in <- file.path(p_dots$dir_plots, paste0('p markers c', data_raw$clust[1], '.png'))
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "Markers by cluster.png")
    file.copy(path_in, path_out, overwrite = TRUE)

    path_in <- file.path(p_dots$dir_plots, 'p_heatmap_all.png')
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "Cluster heat map - all.png")
    file.copy(path_in, path_out, overwrite = TRUE)

    path_in <- file.path(p_dots$dir_plots, 'p_heatmap_filter.png')
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "Cluster heat map - stable.png")
    file.copy(path_in, path_out, overwrite = TRUE)

    #path_in <- file.path(p_dots$dir_plots, paste0('tsne - cluster ',data_raw$clust[1], '.png'))
    #path_out <- file.path(dir_proj, 'exp', "Cluster on t-SNE.png")
    #file.copy(path_in, path_out, overwrite = TRUE)

    path_in <- file.path(str_replace(p_dots$dir_plots, paste0("all_u-", data_raw$stim[1]), "all_u-all"),
                         paste0('umap_clust_grid.png'))
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "UMAP grid.png")
    file.copy(path_in, path_out, overwrite = TRUE)

    path_in <- file.path(p_dots$dir_plots, paste0("p ", p_dots$resp_type, " prog c ", data_raw$clust[1], " fixed.png"))
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "Progressor boxplot - fixed.png")
    file.copy(path_in, path_out, overwrite = TRUE)

    path_in <- file.path(p_dots$dir_plots, paste0("p ", p_dots$resp_type, " prog c ", data_raw$clust[1], " free.png"))
    file.exists(path_in)
    path_out <- file.path(dir_proj, 'exp', "Progressor boxplot - free.png")
    file.copy(path_in, path_out, overwrite = TRUE)
  }

  # =============================
  # Plot - box plots
  # =============================

  # create data for plotting box plots, getting
  # categories for boxplot grouping
  data_plot_box <- data_plot %>%
    mutate(`Time to TB` = ifelse(tfmttb != 0,
                                 paste0("< ", max_ttb, " days"),
                                 Progressor),
           `Time to TB` = ifelse(`Time to TB` == "yes",
                                 paste0("> ", max_ttb, " days"),
                                 `Time to TB`),
           `Time to TB` = ifelse(Progressor == 'no',
                                 'Unknown (control)',
                                 `Time to TB`))

  # create base box plot
  p_box <- ggplot(data_plot_box,
                  aes(x = `Time to TB`, y = resp, col = `Time to TB`)) +
    cowplot::theme_cowplot(font_size = cowplot_theme_font_size)

  # add hor line at y=0
  if(p_dots$assay %in% c('cytof')){
    p_box <- p_box +
      geom_hline(yintercept = 0, linetype = 'dotted', size = 1,
                 col = 'gray50')
  }

  # add points
  if(prop_ind){
    p_box <- p_box + ggforce::geom_sina(aes(size = n_cell),
                                        alpha = 0.8)
  } else p_box <- p_box + ggforce::geom_sina(alpha = 0.8)

  # create break vec so that progressor groups are ordered from close progressors to controls
  break_vec <- unique(data_plot_box$`Time to TB`)
  break_vec <- c( break_vec[str_detect(break_vec, "Unknown")],
                  break_vec[str_detect(break_vec, ">")],
                  break_vec[str_detect(break_vec, "<")])

  p_box <- p_box +
    scale_x_discrete(limits = break_vec)

  # add boxplots
  p_box <- p_box +
    geom_boxplot(alpha = 0.2, outlier.size = -1)

  # adjust colours
  p_box <- p_box +
    scale_colour_manual(values = rev(RColorBrewer::brewer.pal(n = 4, name = 'Oranges')[2:4]),
                        guide = FALSE)

  # add labels
  p_box <- p_box +
    labs(x = "Days to TB", title = "Response by progressor status and days to TB",
         y = y_lab)

  # if points were scaled by size, adjust the legend
  if(prop_ind){
    p_box <- p_box +
      #theme(axis.title.y = element_text(size = y_axis_title_font_size)) +
      scale_size_continuous(name = paste0("Total cell count"),
                            trans = 'log2',
                            range = c(0.1, 3))
  }

  # transform if required
  if('trans_y' %in% names(p_dots)){
    p_box <- p_box +
      scale_y_continuous(name = y_lab,
                         trans = p_dots$trans_y)
  }

  # add grid
  p_box <- p_box +
    cowplot::background_grid(major = 'y', minor = 'y')


  # =============================
  # Plot - progressors over time
  # =============================

  # filter to only choose samples closer to TB
  data_plot_long <- data_plot %>%
    filter(Progressor == 'yes') %>%
    left_join(p_dots$clinical_data %>%
                group_by(SampleID) %>%
                slice(1) %>%
                ungroup() %>%
                select(timeToTBFromVisit, SubjectID, VisitType),
              by = c("SubjectID", "VisitType"))

  # get colour of progressors
  prog_col_ind <- 7
  prog_col <- RColorBrewer::brewer.pal(n = prog_col_ind + 1, name = 'Oranges')[prog_col_ind]

  # create base plot
  p_long <- ggplot(data = data_plot_long,
                   aes(x = timeToTBFromVisit, y = resp))  +
    cowplot::theme_cowplot(font_size = cowplot_theme_font_size)

  # add hor line at y=0
  if(p_dots$assay %in% c('cytof')){
    p_long <- p_long +
      geom_hline(yintercept = 0, linetype = 'dotted', size = 1,
                 col = 'gray50')
  }

  # reverse axis
  p_long <- p_long + scale_x_reverse()

  # add points (with n_cell-dependent size if need be)
  if(prop_ind){
    p_long <- p_long + geom_point(aes(size = n_cell),
                                  col = prog_col,
                                  alpha = 0.8)
  } else p_long <- p_long + geom_point(col = prog_col, alpha = 0.8)

  # add smooth
  p_long <- p_long +
    geom_smooth(formula = y ~ x,
                method = 'lm',
                col = prog_col)

  # if points were scaled by size, adjust the legend
  if(prop_ind){
    p_long <- p_long +
      #theme(axis.title.y = element_text(size = y_axis_title_font_size)) +
      scale_size_continuous(name = paste0("Total cell count"),
                            trans = 'log2',
                            range = c(0.1, 3))
  }

  # transform if required
  if('trans_y' %in% names(p_dots)){
    p_long <- p_long +
      scale_y_continuous(name = y_lab,
                         trans = p_dots$trans_y)
  }

  # add titles
  p_long <- p_long + labs(x = "Days to TB",
                          y = y_lab,
                          title = "Response by days to TB amongst all progressors")

  p_long <- p_long +
    geom_line(data = data_plot_long,
              aes(x = timeToTBFromVisit, y = resp,
                  group = SubjectID),
              alpha = 0.5)

  # add grid
  p_long_prog <- p_long +
    cowplot::background_grid(major = 'xy', minor = 'y')

  # =============================
  # Plot - non-progressors over time
  # =============================

  # filter to only choose samples closer to TB
  data_plot_long <- data_plot %>%
    filter(Progressor == 'no') %>%
    left_join(p_dots$clinical_data %>%
                group_by(SampleID) %>%
                slice(1) %>%
                ungroup() %>%
                select(timeToTBFromVisit, SubjectID, VisitType),
              by = c("SubjectID", "VisitType"))

  # get colour of progressors
  prog_col_ind <- 7
  prog_col <- RColorBrewer::brewer.pal(n = prog_col_ind + 1, name = 'Oranges')[prog_col_ind]

  # create base plot
  p_long <- ggplot(data = data_plot_long,
                   aes(x = DaysSinceEntry, y = resp))  +
    cowplot::theme_cowplot(font_size = cowplot_theme_font_size)

  # add hor line at y=0
  if(p_dots$assay %in% c('cytof')){
    p_long <- p_long +
      geom_hline(yintercept = 0, linetype = 'dotted', size = 1,
                 col = 'gray50')
  }

  # reverse axis
  #p_long <- p_long + scale_x_reverse()

  # add points (with n_cell-dependent size if need be)
  if(prop_ind){
    p_long <- p_long + geom_point(aes(size = n_cell),
                                  col = 'dodgerblue',
                                  alpha = 0.8)
  } else p_long <- p_long + geom_point(col = 'dodgerblue', alpha = 0.8)

  # add smooth
  p_long <- p_long +
    geom_smooth(formula = y ~ x,
                method = 'lm',
                col = 'dodgerblue')

  # if points were scaled by size, adjust the legend
  if(prop_ind){
    p_long <- p_long +
      #theme(axis.title.y = element_text(size = y_axis_title_font_size)) +
      scale_size_continuous(name = paste0("Total cell count"),
                            trans = 'log2',
                            range = c(0.1, 3))
  }

  # transform if required
  if('trans_y' %in% names(p_dots)){
    p_long <- p_long +
      scale_y_continuous(name = y_lab,
                         trans = p_dots$trans_y)
  }

  # add titles
  p_long <- p_long + labs(x = "Days since study entry",
                          y = y_lab,
                          title = "Response by days since study entry amongst all non-progressors")

  # add individual trajectory lines
  p_long <- p_long +
    geom_line(data = data_plot_long,
              aes(x = DaysSinceEntry, y = resp,
                  group = SubjectID),
              alpha = 0.5)

  # add grid
  p_long_ctrl <- p_long +
    cowplot::background_grid(major = 'xy', minor = 'y')

  p_summed <- cowplot::plot_grid(p_box, p_long_prog, p_long_ctrl, ncol = 1)

  # ===========================
  # Save plots
  # ===========================

  gamlsspipeline::save_objects(obj_list = list("Summed response" = p_summed),
                               dir_proj = dir_proj,
                               dir_sub = 'exploration',
                               empty = FALSE,
                               gg_device = '.png',
                               width = 9.21 * 2.5, height = 8 * 2.5 * 1.4)

  # ===========================
  # Cytokine combinations - box plot
  # ===========================

  if(p_dots$assay %in% c('cytof') && p_dots$response_type %in% 'count' && !p_dots$phenotype ){

    # data for plot
    # -------------------

    data_plot <- p_dots$stats_combn_tbl %>%
      mutate(dataset = p_dots$dataset_name,
             cyt_response_type_grp = p_dots$cyt_response_type_grp_curr)

    # filter
    if(length(p_dots$filter) != 0){
      for(i in seq_along(p_dots$filter)){
        var <- names(p_dots$filter)[i]
        val <- p_dots$filter[[i]]
        data_plot <- data_plot[data_plot[[var]] %in% val,]
      }
    }


    # rename freq as resp
    data_plot %<>%# rename(resp = freq_bs)
      mutate(resp = (count_stim/n_cell_stim - count_uns/n_cell_uns) * 1e2)

    # calculate cytokine combination label
    # chnl_vec <- c("Ho165Di", "Gd158Di", "Nd146Di", "Dy164Di", 'Nd150Di', 'Gd156Di')
    #chnl_lab_vec <- setNames(c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6"), chnl_vec)
    # for(chnl in chnl_vec){
    #  if(!chnl %in% names(data_plot)) next
    #  data_plot[[chnl]] <- map_chr(data_plot[[chnl]], function(ind){
    #    ifelse(ind == 1, paste0(chnl_lab_vec[chnl], "+"), paste0(chnl_lab_vec[chnl], "-"))
    #  })
    #}

    #cyt_combn_vec <- rep("", nrow(data_plot))
    #for(chnl in chnl_vec){
    #  if(!chnl %in% names(data_plot)) next
    #  cyt_combn_vec <- paste0(cyt_combn_vec, data_plot[[chnl]])
    #}
    #data_plot[['cyt_combn']] <- cyt_combn_vec

    data_plot %<>% mutate(cyt_combn = str_remove_all(cyt_combn, "~"))
    # data_plot %<>%
    #  mutate(cyt_combn = str_replace(cyt_combn, "Ho165Di", "IFNg") %>%
    #          str_replace("Nd146Di", "TNF") %>%
    #          str_replace("Gd158Di", "IL2") %>%
    #          str_replace("Dy164Di", "IL17") %>%
    #          str_replace("Gd156Di", "IL6") %>%
    #          str_replace("Nd150Di", "IL22") %>%
    #           str_remove_all("~"))
    #for(chnl in chnl_vec){
    #  if(chnl %in% names(data_plot)) data_plot <- data_plot[,-which(colnames(data_plot) == chnl)]
    #}

    # Add clinical data info
    data_plot %<>%
      #mutate(SubjectID = str_sub(SampleID, end = 6),
      #       VisitType = str_sub(SampleID, start = 8)) %>%
      rename(n_cell = n_cell_stim) %>%
      #mutate(fcs = str_remove(fcs, "-time_cut")) %>%
      #left_join(cytof_fcs_to_clin_map %>%
      #            mutate(fcs = stringr::str_remove(OrigFCSName, "-singlets_cleaned") %>%
      #                     stringr::str_remove("-2.1GB_")),
      #          by = 'fcs') %>%
      select(SubjectID, VisitType, cyt_combn, n_cell, resp) %>%
      left_join(clinical_data %>%
                  select(SubjectID, Progressor, timeToTBFromVisit, VisitType),
                by = c("SubjectID", "VisitType"))

    # calculate tfmttb
    data_plot %<>%
      mutate(tfmttb = ifelse(Progressor == 'no', # make zero non-progressor
                                          0,
                                          -9999),
             tfmttb = ifelse(.data$Progressor == 'yes', # make 0 if time to TB is greater than max_ttb; otherwise make
                                          # how many days it is from collapse time to the actual time
                                          pmax(0, .env$max_ttb - .data$timeToTBFromVisit),
                                          .data$tfmttb))

    data_plot %<>%
      group_by(cyt_combn) %>%
      mutate(n_high = sum(resp > 0.02)) %>%
      filter(n_high/n() > 0.1) %>%
      ungroup() %>%
      select(-n_high)

    if(nrow(data_plot) > 1){
      n_cyt_combn <- length(unique(data_plot$cyt_combn))

      # data_plot %<>% filter(!is.na(resp))
      # if(prop_ind) data_plot %<>% mutate(resp = resp/n_cell * 100)

      # create data for plotting box plots, getting
      # categories for boxplot grouping
      data_plot_box <- data_plot %>%
        mutate(`Time to TB` = ifelse(tfmttb != 0,
                                     paste0("< ", max_ttb, " days"),
                                     Progressor),
               `Time to TB` = ifelse(`Time to TB` == "yes",
                                     paste0("> ", max_ttb, " days"),
                                     `Time to TB`),
               `Time to TB` = ifelse(Progressor == 'no',
                                     'Unknown (control)',
                                     `Time to TB`))

      # create base box plot
      p_box <- ggplot(data_plot_box,
                      aes(x = `Time to TB`, y = resp, col = `Time to TB`)) +
        cowplot::theme_cowplot(font_size = cowplot_theme_font_size)

      # add hor line at y=0
      if(p_dots$assay %in% c('cytof')){
        p_box <- p_box +
          geom_hline(yintercept = 0, linetype = 'dotted', size = 1,
                     col = 'gray50')
      }

      # add points
      if(p_dots$assay %in% 'cytof'){
        p_box <- p_box + ggforce::geom_sina(aes(size = n_cell),
                                            alpha = 0.4)
      } else p_box <- p_box + ggforce::geom_sina(alpha = 0.8)

      # create break vec so that progressor groups are ordered from close progressors to controls
      break_vec <- unique(data_plot_box$`Time to TB`)
      break_vec <- c( break_vec[str_detect(break_vec, "Unknown")],
                      break_vec[str_detect(break_vec, ">")],
                      break_vec[str_detect(break_vec, "<")])

      p_box <- p_box +
        scale_x_discrete(limits = break_vec)

      # add boxplots
      p_box <- p_box +
        geom_boxplot(alpha = 0.2, outlier.size = -1)

      # adjust colours
      p_box <- p_box +
        scale_colour_manual(values = rev(RColorBrewer::brewer.pal(n = 4, name = 'Oranges')[2:4]))

      # add labels
      p_box <- p_box +
        labs(x = "Days to TB", title = "Response by progressor status and days to TB",
             y = y_lab)

      # if points were scaled by size, adjust the legend
      if(prop_ind){
        p_box <- p_box +
          #theme(axis.title.y = element_text(size = y_axis_title_font_size)) +
          scale_size_continuous(name = paste0("Total cell count"),
                                trans = 'log2',
                                range = c(0.05, 1))
      }

      # add grid
      p_box <- p_box +
        cowplot::background_grid(major = 'y', minor = 'y') +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

      # facet
      p_box_free <- p_box +
        facet_wrap(~cyt_combn, ncol = 3, scales = 'free')

      p_box_fixed <- p_box +
        facet_wrap(~cyt_combn, ncol = 3, scales = 'fixed') +
        scale_y_continuous(limits = function(lims) c(lims[1], min(lims[2], 0.75)))

      # ===========================
      # Save plots
      # ===========================

      #height_base <- 8 * 2.5 * 1.5
      #height_mult <- 1 + (length(str_split(dataset_name, "_")[[1]]) - 2) * (2.2 - 1)/3
      #height <- height_base * height_mult
      height <- ceiling(n_cyt_combn/3 * 4)
      gamlsspipeline::save_objects(obj_list = list("Cytokine combinations - fixed y-axis scale" = p_box_fixed,
                                                   "Cytokine combinations - free y-axis scale" = p_box_free),
                                   dir_proj = dir_proj,
                                   dir_sub = 'exploration',
                                   empty = FALSE,
                                   gg_device = '.png',
                                   width = 9.21 * 2.5, height = height)
    }



  }
  # ===========================
  # COMPASS
  # ===========================

  if(p_dots$response_type == 'compass'){


    data_meta <- p_dots$compass_obj$data$meta
    #data_meta_orig <- data_meta
    data_meta_add <- p_dots$clinical_data %>%
      filter(SampleID %in% paste0(data_meta$SampleID)) %>%
      group_by(SampleID) %>%
      slice(1) %>%
      ungroup() %>%
      select(SubjectID, VisitType, Progressor, timeToTBFromVisit)

    data_meta_add %<>%
      mutate(ttb = ifelse(is.na(timeToTBFromVisit),
                          max(data_meta_add$timeToTBFromVisit,
                              na.rm = TRUE),
                          timeToTBFromVisit))

    data_meta_add %<>%
      mutate(ttb_cat = ifelse(ttb < 180, "000d - 180d", "99"),
             #ttb_cat = ifelse(ttb >= 90 & ttb < 180, "090d - 180d", ttb_cat),
             ttb_cat = ifelse(ttb >= 180 & ttb < 360, "180d - 360d", ttb_cat),
             ttb_cat = ifelse(ttb >= 360 & ttb < 540, "360d - 540d", ttb_cat),
             ttb_cat = ifelse(ttb >= 540, "540d - 720d", ttb_cat)) %>%
      mutate(ttb_cat = ifelse(Progressor == 'no', '>720d', ttb_cat)) %>%
      mutate(prog = ifelse(Progressor == 'yes', "Case", "Control")) %>%
      mutate(prog_ttb_cat = paste0(prog, ": ", ttb_cat)) %>%
      mutate(`Progressor: Time to TB` = prog_ttb_cat) %>%
      select(-Progressor)

    data_meta %<>%
      left_join(data_meta_add,
                by = c("SubjectID", "VisitType"))

    p_dots$compass_obj$data$meta <- data_meta

    fn_p <- stringr::str_remove(p_dots$html_title, ".html")

    png(file.path(dir_proj, 'exploration', paste0(fn_p, ".png")), width = 400, height = 1000)
    try(suppressMessages(plot(p_dots$compass_obj, 'Progressor: Time to TB')))
    dev.off()

    pdf(file.path(dir_proj, 'exploration', paste0(fn_p, ".pdf")), width = 6, height = 10)
    try(suppressMessages(plot(p_dots$compass_obj, 'Progressor: Time to TB')))
    dev.off()
  }

  if('sig6gene_CorScore' %in% p_dots$conf){
    p_box <- ggplot(data_mod, aes(x = Progressor, y = sig6gene_CorScore, col = Progressor)) +
      background_grid(major = 'xy') +
      geom_boxplot() +
      ggforce::geom_sina() +
      scale_colour_manual(values = c("yes" = "orange", "no" = "dodgerblue")) +
      labs(y = "RISK6 score")

    if(p_dots$assay %in% c('cytof', 'cytof-faust', 'cytof_flowsom') && p_dots$response_type == 'count'){
      data_plot_vs_resp <- data_mod %>% mutate(resp = resp/n_cell * 1e2)
    } else data_plot_vs_resp <- data_mod
    p_verse_resp <- ggplot(data_plot_vs_resp, aes(x = resp, y = sig6gene_CorScore, col = Progressor)) +
      background_grid(major = 'xy') +
      geom_point() +
      geom_smooth(method = 'loess', formula = y ~ x) +
      scale_colour_manual(values = c("yes" = "orange", "no" = "dodgerblue")) +
      labs(x = "Response", y = "RISK6 score")

    p_long_prog <- ggplot(data_mod %>%
                            filter(Progressor == 'yes') %>%
                            mutate(ttb = p_dots$max_ttb - tfmttb * ifelse(p_dots$scale, 1e2, 1)),
                          aes(x = ttb, y = sig6gene_CorScore, col = Progressor)) +
      background_grid(major = 'xy') +
      geom_point() +
      geom_smooth(method = 'loess', formula = y ~ x) +
      scale_colour_manual(values = c("yes" = "orange", "no" = "dodgerblue")) +
      labs(y = "RISK6 score", x = "Days to TB") +
      scale_x_reverse()

    p_long_ctrl <- ggplot(data_mod %>%
                            filter(Progressor == 'no'),
                          aes(x = DaysSinceEntry, y = sig6gene_CorScore, col = Progressor)) +
      background_grid(major = 'xy') +
      geom_point() +
      geom_smooth(method = 'lm', formula = y ~ ns(x)) +
      scale_colour_manual(values = c("yes" = "orange", "no" = "dodgerblue")) +
      labs(y = "RISK6 score", x = "Days since entry")


    gamlsspipeline::save_objects(obj_list = list("Risk6 by prog" = p_box,
                                                 "Risk6 vs resp" = p_verse_resp,
                                                 "Risk6 long prog" = p_long_prog,
                                                 "Risk6 long_ctrl" = p_long_ctrl),
                                 dir_proj = dir_proj,
                                 dir_sub = 'exploration',
                                 empty = FALSE,
                                 gg_device = '.png',
                                 width = 9.21 * 2.5, height = 8)

  }

  # save cluster stability, if available
  if('stability_clust' %in% names(p_dots)){
    saveRDS(p_dots$stability_clust, file.path(dir_proj, 'exploration', 'stability.rds'))
  }


  invisible(TRUE)
}
