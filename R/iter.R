#' @export
#'
#' @param remove_non_tfmttb_int_n logical. If \code{TRUE}, then
#' \code{var_int} cannot be \code{NULL} when the
#' spline explanatory variable is not \code{tfmttb}. Default
#' is \code{TRUE}.
get_iter_tbl <- function(iter_list, remove_non_tfmttb_int_n = TRUE) {

  # cross factors
  iter_tbl <- iter_list %>%
    datautils::cross_df_safe()

  # have only one ttb_max
  # level for non-tfmttb models
  # and set it equal to "-"
  # ------------------
  iter_tbl <- iter_tbl %>%
    dplyr::mutate(
      var_exp_spline_nm = purrr::map(
        var_exp_spline,
        names
      ),
      is_ttb_model = purrr::map_lgl(
        var_exp_spline_nm,
        function(x) "tfmttb" %in% x
      )
    )

  # make ttb_max only have one value if it's not tfmttb as the
  # explanatory variable
  non_ttb_cn_vec <- colnames(iter_tbl) %>%
    setdiff(c("ttb_max"))

  iter_tbl <- iter_tbl %>%
    dplyr::filter(!is_ttb_model) %>%
    dplyr::group_by_at(non_ttb_cn_vec) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(ttb_max = NA_real_) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(
      iter_tbl %>%
        dplyr::filter(is_ttb_model)
    )
  iter_tbl <- iter_tbl %>%
    dplyr::select(-var_exp_spline_nm)

  # interaction
  # ------------------

  if ("var_int" %in% names(iter_list)) {
    # remove it for tfmttb and prog model
    iter_tbl <- iter_tbl %>%
      dplyr::filter(!(is_ttb_model & var_int))

    # set it equal to var_exp and names of var_exp_spline
    var_int_list <- purrr::map(seq_len(nrow(iter_tbl)), function(i) {
      var_int <- iter_tbl$var_int[i]
      if (!var_int) {
        return(NULL)
      }
      c(
        iter_tbl$var_exp[i],
        gsub("^tc~\\w+~", "", names(iter_tbl$var_exp_spline[[i]]))
      )
    })

    # add it to tbl
    iter_tbl <- iter_tbl %>%
      dplyr::mutate(var_int = var_int_list)

    # force interaction terms-only
    # when exp_s is not tfmttb
    if (remove_non_tfmttb_int_n) {
      iter_tbl <- iter_tbl %>%
        dplyr::filter(!(purrr::map_lgl(var_int, is.null) &
          !is_ttb_model))
    }
  }

  # only select one combination
  # for ifng_freq
  if (FALSE) {
    iter_tbl <- iter_tbl %>%
      dplyr::filter(
        !(purrr::map_lgl(
          names(var_exp_spline),
          function(nm) {
            identical(as.character(nm), "tc~flow_ifng~cd4_ifng_freq")
          }
        ) &
          (ds != "cd4_th1_il17" |
            cyt_response_type != "summed" |
            purrr::map_lgl(var_int, function(x) !is.null(x)) |
            stim != "mtb" |
            purrr::map_lgl(var_conf, function(x) !x[[1]] == "none") |
            wins))
      ) %>%
      dplyr::mutate(
        var_exp = ifelse(purrr::map_lgl(
          names(var_exp_spline),
          function(nm) {
            identical(as.character(nm), "tc~flow_ifng~cd4_ifng_freq")
          }
        ),
        "cd4_ifng_freq",
        var_exp
        ),
        var_exp_spline = ifelse(purrr::map_lgl(
          names(var_exp_spline),
          function(nm) {
            identical(as.character(nm), "tc~flow_ifng~cd4_ifng_freq")
          }
        ),
        "none",
        var_exp_spline
        )
      )
  }

  # add boundary knots to tfmttb combinations
  iter_tbl <- set_boundary_knots(iter_tbl)

  # remove helping columns
  iter_tbl <- iter_tbl[, !colnames(iter_tbl) %in% c(
    "is_ttb_model",
    "var_exp_spline_nm"
  )]

  iter_tbl
}


# fix boundary knots if ttb_max is 630
set_boundary_knots <- function(iter_tbl) {
  ttb_model_vec_lgl <- iter_tbl$is_ttb_model
  if (!any(ttb_model_vec_lgl)) {
    return(iter_tbl)
  }
  ttb_model_vec_ind <- which(ttb_model_vec_lgl)
  two_inner_knots_vec_lgl <- purrr::map_lgl(
    ttb_model_vec_ind,
    function(i) {
      iter_tbl$var_exp_spline[[i]]$tfmttb$params$knots %>%
        length() %>%
        magrittr::equals(2)
    }
  )
  if (!any(two_inner_knots_vec_lgl)) {
    return(iter_tbl)
  }
  set_bk_vec_ind <- ttb_model_vec_ind[two_inner_knots_vec_lgl]
  var_exp_spline_list <- iter_tbl %>%
    dplyr::pull(var_exp_spline)
  for (i in set_bk_vec_ind) {
    if (iter_tbl$ttb_max[i] == 630) {
      rep_list <- var_exp_spline_list[[i]]
      rep_list$tfmttb$params$Boundary.knots <- c(0.3, 4.5)
      var_exp_spline_list[[i]] <- rep_list
    }
  }
  iter_tbl %>%
    dplyr::mutate(var_exp_spline = var_exp_spline_list)
}

#' @title Return NULL if object is "none" or NULL
#' @export
set_none_to_null <- function(x, each_elem = TRUE) {
  switch(as.character(each_elem),
    "TRUE" = purrr::map(x, .set_none_to_null) %>%
      setNames(names(x)),
    .set_none_to_null(x)
  )
}

.set_none_to_null <- function(x) {
  if (class(try(x[[1]])) == "try-error") {
    return(x)
  }
  switch(as.character(identical(x[[1]], "none") || is.null(x[[1]])),
    "TRUE" = NULL,
    x
  )
}

#' @title Print iterator
#'
#' @description Prints all factors
#' in current iterator that are different to the last.
#' Assumes that \code{iter} is the name of the iterator
#'
#' @param ind integer. Row index of the iterator.
#' @export
print_iter <- function(iter, ind) {
  if (ind == 1) {
    iter_old <- iter

    iter_old[, seq_len(ncol(iter_old))] <- "&$#("
  } else {
    iter_old <- parent.frame(2)$iter_old
    if (!identical(
      colnames(iter),
      colnames(iter_old)
    )
    ) {
      iter_old <- iter
      iter_old[, seq_len(ncol(iter_old))] <- "&$#("
    }
  }
  purrr::walk(seq_len(ncol(iter)), function(j) {
    if (!identical(iter[[j]], iter_old[[j]])) {
      cat("--- ", paste0(names(iter[[j]]), " ---\n"))
      cat(paste0(
        paste0(iter[[j]], collapse = "; "), "\n"
      ))
    }
  })
  iter_old <- iter
  assign("iter_old",
    value = iter_old,
    envir = parent.frame(2)
  )
  invisible(iter_old)
}

#' @export
order_iter <- function(iter, factor_last = NULL) {
  if (is.null(factor_last)) {
    return(iter)
  }
  order_vec_cn <- c(
    setdiff(colnames(iter), factor_last),
    factor_last[factor_last %in% colnames(iter)]
  )
  iter[, order_vec_cn]
}
