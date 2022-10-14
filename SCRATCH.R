winsorise <- function(data_raw, wins, p_dots, iter) {
  if (wins == "wins_n") {
    return(data_raw)
  }

  wins_var <- stringr::str_sub(wins, start = 6)
  wins_var_vec <- purrr::map_chr(seq_len(stringr::str_length(wins_var)), function(i) {
    stringr::str_sub(wins_var, i, i)
  })
  if ("y" %in% wins_var_vec) {
    resp_vec <- data_raw[["resp"]]
    if (iter$var_offset == "n_cell") {
      browser()
    } else if (iter$var_offset != "none") {
      stop(paste0(
        "value for var_offset of ",
        iter$var_offset,
        " not recognised"
      ))
    }
    sd_var <- sd(resp_vec)
    mad_var <- mad(resp_vec)
    max_var <- max(
      mean(resp_vec) + 2 * sd_var,
      median(resp_vec) + 2 * mad_var
    )
    min_var <- min(
      mean(resp_vec) - 2 * sd_var,
      median(resp_vec) - 2 * mad_var
    )
    data_raw[, "resp"] <- pmax(pmin(resp_vec, max_var), min_var)
  }
  if ("x" %in% wins_var_vec) {
    var_vec <- c(iter$var_exp, names(iter$var_exp_s))
    for (var in var_vec) {
      if (!is.numeric(data_raw[[var]])) next
      sd_var <- sd(data_raw[[var]])
      mad_var <- mad(data_raw[[var]])
      max_var <- max(
        mean(data_raw[[var]]) + 3 * sd_var,
        median(data_raw[[var]]) + 3 * mad_var
      )
      min_var <- min(
        mean(data_raw[[var]]) - 3 * sd_var,
        median(data_raw[[var]]) - 3 * mad_var
      )
      data_raw[, var] <- pmax(pmin(data_raw[[var]], max_var), min_var)
    }
  }
  data_raw
}
