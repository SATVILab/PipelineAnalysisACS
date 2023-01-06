#' @title Plot model est and ci for int between cat and num exp var
plot_disp_int_cat_num <- function(mod, .data, data_nm,
                                  var_cat, var_num,
                                  var_offset,
                                  var_dep,
                                  max_sd = NULL,
                                  axis_x_reverse = FALSE,
                                  cat_to_col = NULL,
                                  axis_lab = NULL,
                                  add_test = "lr",
                                  trans = scales::identity_trans(),
                                  dir_test,
                                  table_coord = c(0.1, 1),
                                  table_size_skip = 0.05,
                                  table_size_text = 10,
                                  range_extend = 0,
                                  point_alpha = 1,
                                  limits_include = NULL,
                                  gg_theme = cowplot::theme_cowplot(),
                                  grid = "xy",
                                  point_size = NULL,
                                  txt_parse = FALSE) {
  # prep
  # --------------

  # get estimates and ci's
  eff_obj <- run_effects_Effect(
    data = .data,
    nm = data_nm,
    var = c(var_cat, var_num),
    mod = mod
  )

  plot_tbl_eff <- eff_obj |>
    as.data.frame() |>
    tibble::as_tibble()

  # check if n_cell was used as number of offsets
  n_cell_ind <- switch(as.character(is.null(var_offset)),
    "TRUE" = FALSE,
    "FALSE" = grepl("n_cell", var_offset)
  )

  # make response a frequency if offset
  # is number of cells
  if (n_cell_ind) {
    n_cell <- exp(eff_obj$offset)
    plot_tbl_eff <- plot_tbl_eff |>
      dplyr::mutate(
        fit = fit / n_cell * 1e2,
        lower = lower / n_cell * 1e2,
        upper = upper / n_cell * 1e2
      )
  }

  # prep data_raw
  # -------------------

  plot_tbl_raw <- .data

  # account for offset, if need be
  if (!is.null(var_offset)) {
    plot_tbl_raw[, var_dep] <- plot_tbl_raw[, var_dep] /
      plot_tbl_raw[, var_offset]
    if (n_cell_ind) plot_tbl_raw[, var_dep] <- plot_tbl_raw[, var_dep] * 1e2
  }

  # restrict plot_tbl_eff to where
  # we actually have data
  # ----------------------

  if (!(var_cat == "Progressor" && var_num == "tfmttb")) {
    cat_vec <- unique(plot_tbl_eff[[var_cat]]) |>
      as.character()

    plot_tbl_eff <- purrr::map_df(cat_vec, function(cat) {
      eff_tbl_filter <- plot_tbl_eff |>
        dplyr::filter(.data[[var_cat]] == .env$cat)
      raw_tbl_filter <- plot_tbl_raw |>
        dplyr::filter(.data[[var_cat]] == .env$cat)
      pad_var_num <- diff(range(raw_tbl_filter[[var_num]],
        na.rm = TRUE
      )) * 0.01
      lims_tbl_raw <- raw_tbl_filter |>
        dplyr::summarise(
          min = min(.data[[var_num]], na.rm = TRUE) -
            pad_var_num,
          max = max(.data[[var_num]], na.rm = TRUE) +
            pad_var_num
        )

      eff_tbl_filter_max <- eff_tbl_filter |>
        dplyr::filter(.data[[var_num]] >= lims_tbl_raw$max)
      if (nrow(eff_tbl_filter_max) > 0) {
        max_num_var <- eff_tbl_filter_max |>
          dplyr::filter(.data[[var_num]] == min(.data[[var_num]],
            na.rm = TRUE
          )) |>
          dplyr::pull(var_num)
        max_num_var <- min(max_num_var)
      } else {
        max_num_var <- max(eff_tbl_filter[[var_num]])
      }

      eff_tbl_filter_min <- eff_tbl_filter |>
        dplyr::filter(.data[[var_num]] <= lims_tbl_raw$min)
      if (nrow(eff_tbl_filter_min) > 0) {
        min_num_var <- eff_tbl_filter_min |>
          dplyr::filter(.data[[var_num]] == max(.data[[var_num]],
            na.rm = TRUE
          )) |>
          dplyr::pull(var_num)
        min_num_var <- max(min_num_var)
      } else {
        min_num_var <- min(eff_tbl_filter[[var_num]])
      }
      eff_tbl_filter |>
        dplyr::filter(
          .data[[var_num]] >= min_num_var,
          .data[[var_num]] <= max_num_var
        )
    })
  }


  # recalibrate tfmttb
  # ------------------

  if (var_cat == "Progressor" && var_num == "tfmttb") {
    ttb_max <- rlang::caller_env()$iter$ttb_max

    plot_tbl_raw <- plot_tbl_raw |>
      dplyr::mutate(tfmttb = ttb_max - (tfmttb * 1e2))
  }

  # stretch non-progressor estimates
  # across range of numeric variable
  if (var_cat == "Progressor" && var_num == "tfmttb") {
    eff_tbl_prog <- plot_tbl_eff |>
      dplyr::filter(Progressor == "yes")

    eff_tbl_ctrl <- purrr::map_df(eff_tbl_prog$tfmttb, function(x) {
      plot_tbl_eff |>
        dplyr::filter(
          Progressor == "no",
          tfmttb == 0
        ) |>
        dplyr::mutate(tfmttb = .env$x)
    })

    plot_tbl_eff <- eff_tbl_prog |>
      dplyr::bind_rows(eff_tbl_ctrl)

    ttb_max <- rlang::caller_env()$iter$ttb_max

    plot_tbl_eff <- plot_tbl_eff |>
      dplyr::mutate(tfmttb = ttb_max - (tfmttb * 1e2))
  }

  p <- ggplot(plot_tbl_eff, aes(
    x = .data[[var_num]], y = fit,
    col = .data[[var_cat]],
    fill = .data[[var_cat]]
  )) +
    gg_theme +
    cowplot::background_grid(major = "xy", minor = "y") +
    # scale_x_reverse() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
      alpha = 0.5
    ) +
    geom_line(aes(y = fit))

  cat_to_col <- switch(as.character(is.null(cat_to_col)),
    "TRUE" = suppressWarnings(
      RColorBrewer::brewer.pal(
        n = length(unique(.data[[var_cat]])),
        name = "Set1"
      )
    ),
    "FALSE" = cat_to_col
  )


  p <- p +
    scale_colour_manual(
      values = cat_to_col
    ) +
    scale_fill_manual(
      values = cat_to_col
    )

  if (axis_x_reverse) p <- p + ggplot2::scale_x_reverse()

  if (!is.null(axis_lab)) {
    p <- p +
      labs(x = axis_lab[1], y = axis_lab[2])
  }

  # make raw data plot
  # -------------------

  # format raw data
  point_geom_list <- list(
    data = plot_tbl_raw,
    aes(y = resp),
    alpha = point_alpha
  )
  point_geom_list <- point_geom_list |>
    append(
      switch(as.character(is.null(point_size)),
        "TRUE" = list(),
        "FALSE" = list(size = point_size)
      )
    )

  point_geom <- do.call(
    geom_point,
    point_geom_list
  )

  p_raw <- p + point_geom

  # restrict to between 0 and 100 output if modelling frequencies
  # ------------------
  if (n_cell_ind) {
    max_pt <- max(plot_tbl_raw$resp)
    max_est <- quantile(plot_tbl_eff$upper, 0.75)
    max_val <- max(max_pt, max_est)

    if (!is.null(max_sd)) {
      sd_quant_vec <- c(
        quantile(plot_tbl_raw$resp, 0.025),
        quantile(plot_tbl_raw$resp, 0.975)
      )
      plot_tbl_raw_resp_restr <- plot_tbl_raw$resp[
        plot_tbl_raw$resp >= sd_quant_vec[1] &
          plot_tbl_raw$resp <= sd_quant_vec[2]
      ]
      sd_pt <- sd(plot_tbl_raw_resp_restr)
      mean_pt <- mean(plot_tbl_raw_resp_restr)
      sd_min_max_vec <- c(
        mean_pt - max_sd * sd_pt,
        mean_pt + max_sd * sd_pt
      )
      max_val <- min(max_val, sd_min_max_vec[2])
    }

    p <- p + ggplot2::coord_cartesian(ylim = c(0, max_val))
    p_raw <- p_raw +
      coord_cartesian(ylim = c(0, max_val))
  } else {
    max_pt <- max(plot_tbl_raw$resp)
    max_est <- quantile(plot_tbl_eff$upper, 0.75)
    max_val <- max(max_pt, max_est)

    min_pt <- min(plot_tbl_raw$resp)
    min_est <- quantile(plot_tbl_eff$lower, 0.25)
    min_val <- min(min_pt, min_est)

    if (!is.null(max_sd)) {
      sd_quant_vec <- c(
        quantile(plot_tbl_raw$resp, 0.025),
        quantile(plot_tbl_raw$resp, 0.975)
      )
      plot_tbl_raw_resp_restr <- plot_tbl_raw$resp[
        plot_tbl_raw$resp >= sd_quant_vec[1] &
          plot_tbl_raw$resp <= sd_quant_vec[2]
      ]
      sd_pt <- sd(plot_tbl_raw_resp_restr)
      mean_pt <- mean(plot_tbl_raw_resp_restr)
      sd_min_max_vec <- c(
        mean_pt - max_sd * sd_pt,
        mean_pt + max_sd * sd_pt
      )
      min_val <- max(min_val, sd_min_max_vec[1])
      max_val <- min(max_val, sd_min_max_vec[2])
    }

    p <- p + ggplot2::coord_cartesian(ylim = c(min_val, max_val))
    p_raw <- p_raw +
      coord_cartesian(ylim = c(min_val, max_val))
  }



  p_list <- list(p, p_raw)

  # add a results table
  # ----------------

  if (!is.null(add_test)) {
    test_tbl <- readRDS(file.path(dir_test, paste0(add_test, ".rds")))

    test_tbl <- test_tbl |>
      dplyr::mutate(var = purrr::map_chr(var, function(x) {
        x <- x |>
          stringr::str_replace("Progressor; tfmttb", "Progression") |>
          stringr::str_replace("^Progressor$", "Distal TB") |>
          stringr::str_replace("; ", " and/or ") |>
          stringr::str_replace(":", " int. with ") |>
          stringr::str_replace("^tfmttb$", "Proximal TB") |>
          stringr::str_replace("risk6", "RISK6") |>
          stringr::str_replace("^Progressor$", "Distal TB") |>
          stringr::str_replace("^Progressor ", "Progression ")
      })) |>
      dplyr::mutate(p = signif(p, digits = 3))

    if (!test_tbl$var[1] == "Progression") {
      test_tbl$var[which(test_tbl$var == "Distal TB")] <- "Progression"
    }

    res_text <- paste0(test_tbl$var, ": ", signif(test_tbl$p, digits = 3))

    p_list <- purrr::map(p_list, function(p) {
      UtilsGGSV::add_text_column(
        p = p,
        x = range(c(plot_tbl_eff[[var_num]], plot_tbl_raw[[var_num]])),
        y = switch(as.character(n_cell_ind),
          "TRUE" = c(0, max_val),
          "FALSE" = c(min_val, max_val)
        ),
        font_size = table_size_text,
        coord = c(ifelse(var_num == "tfmttb", 0.95, 0.05), 0.95),
        skip = table_size_skip,
        text = res_text,
        parse = txt_parse
      )
    })
  }

  p_list |>
    setNames(c(
      "p_fit",
      "p_fit_raw"
    ))
}
