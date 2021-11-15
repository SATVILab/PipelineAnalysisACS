.get_np_results <- function(data_mod, p_dots, iter) {
  resp_vec <- data_mod$resp
  if (!is.null(iter$var_offset)) {
    resp_vec <- resp_vec / data_mod[[iter$var_offset]]
  }

  if (names(iter$var_exp_spline) == "tfmttb") {
    results_tbl_np <- tibble::tibble(var = character(0), p = numeric(0))
    var_exp_s <- iter$var_exp_spline
    params_list <- iter[["var_exp_spline"]][["tfmttb"]]$params

    exp_s_vec <- data_mod$tfmttb
    cat_vec_tfmttb <- cut(exp_s_vec,
      breaks = c(
        min(exp_s_vec) - 1,
        params_list$knots,
        max(exp_s_vec) + 1
      )
    )

    p_tfmttb <- kruskal.test(
      x = resp_vec,
      g = cat_vec_tfmttb
    )$p.value
    results_tbl_np <- results_tbl_np %>%
      dplyr::bind_rows(
        tibble::tibble(
          var = "tfmttb",
          p = p_tfmttb
        )
      )

    p_prog <- kruskal.test(
      x = resp_vec,
      g = data_mod$Progressor
    )$p.value
    results_tbl_np <- results_tbl_np %>%
      dplyr::bind_rows(
        tibble::tibble(
          var = "Progressor",
          p = p_prog
        )
      )

    cat_vec_int <- purrr::map_chr(seq_len(nrow(data_mod)), function(i) {
      if (data_mod$Progressor[[i]] == "no") {
        return("no")
      }
      cat_vec_tfmttb[i]
    })
    p_prog_tfmttb <- kruskal.test(
      x = resp_vec,
      g = cat_vec_int
    )$p.value

    results_tbl_np <- results_tbl_np %>%
      dplyr::bind_rows(
        tibble::tibble(
          var = "Progressor; tfmttb",
          p = p_prog_tfmttb
        )
      )

    results_tbl_np <- results_tbl_np[c(3, 2, 1), ]
  } else {
    results_tbl_np <- tibble::tibble(var = character(0), p = numeric(0))
    var_exp_s <- iter$var_exp_spline

    if (!is.null(iter$var_int)) {
      cat_list <- purrr::map(iter$var_int[[1]], function(x) {
        exp_vec <- data_mod[[x]]
        if (x %in% names(iter$var_exp_spline)) {
          params_list <- iter$var_exp_spline[[x]]$params
          if ("df" %in% names(params_list)) {
            df <- params_list$df
            cat_vec <- cut(exp_vec, breaks = df + 1)
          } else if ("knots" %in% names(params_list)) {
            cat_vec <- cut(exp_vec,
              breaks = c(
                min(exp_vec) - 1,
                params_list$knots,
                max(exp_vec) + 1
              )
            )
            return(cat_vec)
          } else {
            stop("neither knots nor df specified for splines term")
          }
          return(cat_vec)
        }
        if (is.numeric(exp_vec)) {
          exp_vec <- cut(exp_vec, breaks = df + 1)
        }
        exp_vec
      })
      cat_vec <- paste0(
        cat_list[[1]],
        cat_list[[2]]
      )
      p <- kruskal.test(
        x = resp_vec,
        g = cat_vec
      )$p.value
      results_tbl_np <- results_tbl_np %>%
        dplyr::bind_rows(
          tibble::tibble(
            var = paste0(iter$var_int, collapse = "; "),
            p = p
          )
        )
    }

    if (!is.null(iter$var_exp)) {
      exp_vec <- data_mod[[iter$var_exp]]
      if (is.numeric(exp_vec)) {
        exp_vec <- cut(exp_vec, breaks = df + 1)
      }

      p <- kruskal.test(
        x = resp_vec,
        g = exp_vec
      )$p.value

      results_tbl_np <- results_tbl_np %>%
        dplyr::bind_rows(
          tibble::tibble(
            var = iter$var_exp,
            p = p
          )
        )
    }

    if (!is.null(iter[["var_exp_spline"]])) {
      df <- var_exp_s[[1]]$params$df
      exp_s_vec <- data_mod[[names(iter$var_exp_spline)]]
      cat_vec <- cut(exp_s_vec, breaks = df + 1)
      p <- kruskal.test(
        x = resp_vec,
        g = cat_vec
      )$p.value
      results_tbl_np <- results_tbl_np %>%
        dplyr::bind_rows(
          tibble::tibble(
            var = names(var_exp_s),
            p = p
          )
        )
    }
  }


  results_tbl_np
}
