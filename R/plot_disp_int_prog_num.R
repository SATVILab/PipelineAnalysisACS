#' @title Plot model est and ci for int between cat and num exp var
plot_disp_int_cat_num <- function(mod, .data, data_nm,
                                  var_cat, var_num,
                                  var_offset,
                                  var_dep,
                                  axis_x_reverse = FALSE,
                                  cat_to_col = NULL, axis_lab = NULL,
                                  add_test = "lr",
                                  dir_test) {

  # prep
  # --------------

  # get non-progressor variable name
  # var_p_n <- setdiff(var_int, "Progressor")
  #var_alt_ind <- which(colnames(eff_tbl) == np_var)
  #colnames(eff_tbl)[var_alt_ind] <- 'var_alt'


  # get estimates and ci's
  eff_obj <- run_effects_Effect(
    data = .data,
    nm = data_nm,
    var = c(var_cat, var_num),
    mod = mod
  )

  eff_tbl <- eff_obj %>%
    as.data.frame() %>%
    tibble::as_tibble()


  # check if n_cell was used as number of offsets
  n_cell_ind <- switch(
    as.character(is.null(var_offset)),
    "TRUE" = FALSE,
    "FALSE" = grepl("n_cell", var_offset))

  # make response a frequency if offset
  # is number of cells
  if(n_cell_ind){
    n_cell <- exp(eff_obj$offset)
    eff_tbl <- eff_tbl %>%
      dplyr::mutate(fit = fit / n_cell * 1e2,
                    lower = lower / n_cell * 1e2,
                    upper = upper / n_cell * 1e2)
  }

  # calculate plot bounds
  # -----------------------

  # stretch non-progressor estimates
  # across range of numeric variable
  if (var_cat == "Progressor" && var_num == "tfmttb") {

  }

  p <- ggplot(eff_tbl, aes(x = .data[[var_num]], y = fit,
                      col = .data[[var_cat]],
                      fill = .data[[var_cat]])) +
    cowplot::background_grid(major = 'xy', minor = 'y') +
    #scale_x_reverse() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5) +
    geom_line(aes(y = fit))

  cat_to_col <- switch(
    as.character(is.null(cat_to_col)),
    "TRUE" = suppressWarnings(
      RColorBrewer::brewer.pal(n = length(unique(.data[[var_cat]])),
                               name = "Set1")
      ),
    "FALSE" = cat_to_col)

  p <- p +
    scale_colour_manual(
      values = cat_to_col
    ) +
    scale_fill_manual(
      values = cat_to_col
    )

  if (axis_x_reverse) p <- p + scale_x_reverse()

  if (!is.null(axis_lab)) {
    p <- p +
      labs(x = axis_lab[1], y = axis_lab[2])
  }

  # make raw data plot
  # -------------------

  # format raw data

  plot_tbl_raw <- .data

  # account for offset, if need be
  if (!is.null(var_offset)) {
    plot_tbl_raw[, var_dep] <- plot_tbl_raw[, var_dep] / plot_tbl_raw[, var_offset]
    if (n_cell_ind) plot_tbl_raw[, var_dep] <- plot_tbl_raw[, var_dep] * 1e2
  }

  p_raw <- p +
    geom_point(data = plot_tbl_raw,
               aes(y = resp))

  # restrict to between 0 and 100 output if modelling frequencies
  # ------------------
  if (n_cell_ind) {
    max_pt <- max(plot_tbl_raw$resp) + 0.01 * diff(range(plot_tbl_raw$resp))
    max_est <- quantile(eff_tbl$upper, 0.75)
    max_val <- min(100, max_pt, max_est)
    p <- p + coord_cartesian(ylim = c(0, max_val))
    p_raw <- p_raw +
      coord_cartesian(ylim = c(0, max_val))
  }

  p_list <- list(p, p_raw)

  # add a results table
  # ----------------

  if(!is.null(add_test)) {
    test_tbl <- readRDS(file.path(dir_test, paste0(add_test, ".rds")))

    test_tbl <- test_tbl %>%
      dplyr::mutate(var = purrr::map_chr(var, function(x){
        x %>%
          stringr::str_replace("; ", " and/or ") %>%
          stringr::str_replace(":",  " int. with ")
      }))

    res_text <- paste0(test_tbl$var, ": ", signif(test_tbl$p, digits = 3))

    p_list <- purrr::map(p_list, function(p){
      p <- cowplot::ggdraw(p)
      for(i in seq_along(res_text)){
        p <- p +
          cowplot::draw_text(res_text[i], x = 0.1675, y = 0.95 - 0.05 * (i-1),
                             size = 12 , hjust = 0)
      }
      p
    })
  }

  p_list %>%
    setNames(c("p_fit_manual_without_raw_data",
               "p_fit_manual_with_raw_data"))
}
