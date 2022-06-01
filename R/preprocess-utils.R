fortify_vars <- function(iter,
                         regex = "var_*",
                         cn = NULL,
                         default_ds = "data_raw") {
  stopifnot(is.list(iter))
  if (!is.null(regex)) {
    iter_nm_vec <- switch(as.character(is.data.frame(iter)),
      "TRUE" = colnames(iter),
      "FALSE" = names(iter)
    )
    cn <- c(
      cn,
      iter_nm_vec[grepl(regex, iter_nm_vec)]
    )
  }
  cn <- unique(cn)
  if (length(cn) == 0) {
    return(iter)
  }

  for (.cn in cn) {
    print(.cn)
    iter[[.cn]] <- .fortify_var(
      elem = iter[[.cn]],
      default_ds = default_ds
    )
  }
  iter
}

.fortify_var <- function(elem, default_ds) {
  if (is.null(elem)) return(list(NULL))
  if (identical(list(NULL), elem)) return(elem)
  if (!is.list(elem) && is.vector(elem)) {
    elem_out <- lapply(elem, function(x) {
      list(
        nm = x,
        cn_orig = x,
        cn_new = x,
        ds = default_ds,
        trans = NULL,
        spline = NULL
      )
    })
  } else if (is.list(elem)) {
    if (any(vapply(elem, Negate(is.list), FUN.VALUE = logical(1)))) {
      stop(
        "Variables must be specified as a list of lists",
        call. = FALSE
      )
    }
    elem_out <- lapply(elem, function(x) {
      if (!("cn_orig" %in% names(x))) {
        stop("cn_orig not supplied")
      }
      if (!"nm" %in% names(x)) {
        x[["nm"]] <- setNames(unlist(x[["cn_orig"]]), NULL)
      }
      if (!"cn_new" %in% names(x)) {
        x[["cn_new"]] <- setNames(unlist(x[["cn_orig"]]), NULL)
      }
      if (!"ds" %in% names(x)) x[["ds"]] <- default_ds
      if (!"trans" %in% names(x)) x <- append(x, list("trans" = NULL))
      if (!"spline" %in% names(x)) x <- append(x, list("spline" = NULL))
      if (!is.null(x[["spline"]])) {
        stopifnot("fn" %in% names(x$spline))
        if ("args" %in% names(x$spline)) stopifnot(is.list(x$spline$args))
        extra_elems <- setdiff(names(x$spline), c("pkg", "fn", "args"))
        if (length(extra_elems) > 0) {
          stop(
            paste0(
              "The following are extra arguments to spline: ",
              paste0(extra_elems, collapse = "; ")
            )
          )
        }
      }
      if (!is.null(x$trans)) stopifnot(is.function(x$trans))
      x[c("nm", "cn_orig", "cn_new", "ds", "trans", "spline")]
    })
  } else {
    stop("elem not either a list or a vector (or NULL)")
  }
  elem_out
}

#' @export
filter_using_list <- function(.data, filter_list) {
  for (i in seq_along(filter_list)) {
    nm <- names(filter_list)[i]
    value <- filter_list[[i]]
    .data <- .data[.data[[nm]] %in% value, ]
  }
  .data
}

get_tfmttb <- function(.data, ttb_max, ttb_min) {

  if (missing(ttb_max)) stop("ttb_max must be specified")
  if (missing(ttb_min)) stop("ttb_min must be specified")

  if (!is.null(ttb_min)) {
    .data <- .data %>%
      dplyr::mutate(
        timeToTBFromVisit = pmax(.data$timeToTBFromVisit, ttb_min)
      )
  }

  ttb_max <- ifelse(
    is.null(ttb_max) || is.na(ttb_max),
    max(.data$timeToTBFromVisit,na.rm = TRUE),
    ttb_max
  )

  .data %>%
    dplyr::mutate(
      tfmttb = ifelse(Progressor == "no", # make zero non-progressor
        0,
        -9999
      ),
      tfmttb = ifelse(.data$Progressor == "yes",
        pmax(0, ttb_max - .data$timeToTBFromVisit),
        .data$tfmttb
      )
    )
}

trans <- function(.data, trans) {
  if (is.null(trans)) {
    return(.data)
  }
  if (is.function(trans)) {
    .data <- .data %>%
      dplyr::mutate(resp = trans(resp))
    return(.data)
  }
  if (!is.character(trans)) {
    stop("trans must be NULL, a function or a character")
  }

  trans_fn <- switch(trans,
    "min_max" = ,
    "within_01" = function(x) {
      range_x <- range(x)
      diff_x <- diff(range_x)
      # now (not strictly) within 0,1
      x <- (x - range_x[1]) / diff_x
      # strictly within 01:
      small_diff <- diff_x * 0.005
      pmax(small_diff, pmin(x, range_x[2] - small_diff))
    },
    "nudge_off_01" = function(x) {
      range_x <- range(x)
      diff_x <- diff(range_x)
      max_x <- max(1 - diff_x, max(x[x < 1]))
      min_x <- min(diff_x, min(x[x > 0]))
      pmax(min_x, pmin(x, max_x))
    },
    "pos" = function(x) pmax(0, x),
    "pos_strict" = function(x) {
      range_x <- diff(range(x))
      diff_x <- range_x * 0.01
      min_x <- min(diff_x, min(x[x > 0]))
      pmax(min_x, x)
    },
    "sqrt" = sqrt,
    "sqrt_pos" = function(x) sqrt(pmax(0, x)),
    stop("trans supplied but not recognised")
  )
  .data %>%
    dplyr::mutate(resp = trans_fn(resp))
}

#' @title Calculate bs-freq for CyTOF phenotypes
#'
#' @param .data dataframe.
#' @param cyt_response_type "summed" or "combn". Specifies
#' type of response to create.
#' Should be able to specify "single" as well.
#' @export
prep_bs_freq <- function(.data, cyt_response_type) {
  .data <- .data %>%
    dplyr::filter(grepl("\\+", .data$cyt_combn))

  # sum over cytokines
  .data <- switch(cyt_response_type,
    "summed" = .data %>%
      UtilsCytoRSV::sum_over_markers(
        grp = c(
          "batch_sh", "SubjectID", "VisitType",
          "stim", "n_cell_stim", "n_cell_uns"
        ),
        cmbn = "cyt_combn",
        levels = c("-", "+"),
        resp = c("count_stim", "count_uns")
      ) %>%
      dplyr::mutate(cyt_combn = "summed"),
    "combn" = .data %>%
      dplyr::select(
        batch_sh, SubjectID, VisitType,
        cyt_combn,
        stim, n_cell_stim, n_cell_uns,
        count_stim, count_uns
      ) %>%
      dplyr::filter(grepl("\\+", cyt_combn)),
    stop(paste0(
      paste0(cyt_response_type, collapse = "/"),
      " not available for cyt_response_type in prep_bs_freq"
    ))
  )



  .data
}

#' @export
get_var <- function(iter,
                    iter_cn = c(
                      "var_dep",
                      "var_offset",
                      "var_re",
                      "var_conf",
                      "var_exp",
                      "var_exp_spline"
                      ),
                    var = NULL) {
  # force exact matches
  if (!tibble::is_tibble(iter)) {
    if (is.data.frame(iter)) {
      iter <- tibble::as_tibble(iter)
    } else {
      stop("iter must be a tibble", .call = FALSE)
    }
  }

  # add var_exp_spline entries
  if ("var_exp_spline" %in% iter_cn) {
    var_vec_exp_spline <- vapply(
      iter$var_exp_spline, function(x) x[[1]]$var,
      character(1)
    )
  } else var_vec_exp_spline <- NULL

  # add other variables than var_exp_spline
  iter_cn_vec_non_exp_spline <- setdiff(iter_cn, "var_exp_spline")
  var_vec_other <- NULL
  for (i in seq_along(iter_cn_vec_non_exp_spline)) {
    var_vec_other <- c(
      var_vec_other,
      iter[[iter_cn_vec_non_exp_spline[i]]] %>%
        unlist() %>%
        setNames(NULL)
    )
  }

  unique(c(var_vec_exp_spline, var_vec_other, var))
}

#' @title Get clinical variables to add to a dataset
#'
#' @param iter one-row dataframe.
#' @param iter_cn character vector.
#' Names of columns in
#' \code{iter} to get clinical variable names from.
#' @param var character vector.
#' Variables to add not specificied by iter.
#' @export
get_var_clin <- function(var,
                          var_clin_possible = c(
                            "DaysSinceEntry",
                            "AgeAtLastBirthDay",
                            "tfmttb"
                          )) {
  intersect(var, var_clin_possible)
}

#' @export
#'
#' @param
#'
#' @details
#' - Join on SubjectID and VisitType
#' - Always add Progressor, tfmttb and timeToTBFromVisit
#' - Returns only rows from .data
#' where sample is before TB diagnosis
add_clinical_data_and_filter <- function(.data,
                                         cols_add = NULL,
                                         cols_join = c("SubjectID", "VisitType"),
                                         cols_add_perm = c(
                                           "Progressor",
                                           "timeToTBFromVisit",
                                           "tfmttb"
                                         ),
                                         ttb_min = NULL,
                                         ttb_max) {
  cols_add <- unique(c(cols_add_perm, cols_add)) %>%
    setdiff("n_cell")
  cols_add <- unique(cols_add)
  cols_join <- unique(cols_join)
  clinical_data_add <- TuberculomicsCompendium::clinical_data
  clinical_data_add <- clinical_data_add[, colnames(clinical_data_add) %in%
    c(cols_join, cols_add)]
  clinical_data_add <- clinical_data_add %>%
    dplyr::mutate(SampleID = paste0(.data$SubjectID, "_", .data$VisitType)) %>%
    dplyr::filter(SampleID %in% unique(
      .env$.data %>%
        dplyr::mutate(SampleID = paste0(.data$SubjectID, "_", .data$VisitType)) %>%
        dplyr::pull("SampleID")
    )) %>%
    dplyr::select(-SampleID)

  clinical_data_add <- clinical_data_add %>%
    dplyr::group_by_at(c(cols_join, setdiff(cols_add, "tfmttb"))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # check that there is only one row per key
  # (key is combn of joining col levels)
  clinical_data_add_check <- clinical_data_add %>%
    dplyr::group_by_at(cols_join) %>%
    dplyr::summarise(cnt = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(g_1 = any(cnt > 1)) %>%
    dplyr::pull(g_1)
  if (clinical_data_add_check) stop("clinical_data has more than one row per key")

  clinical_data_add <- clinical_data_add %>%
    dplyr::filter(timeToTBFromVisit >= 0 | is.na(timeToTBFromVisit))

  clinical_data_add <- clinical_data_add %>%
    dplyr::filter(Progressor %in% c("yes", "no"))

  if (!is.null(ttb_min)) {
    clinical_data_add <- clinical_data_add %>%
      dplyr::mutate(timeToTBFromVisit = pmax(.data$timeToTBFromVisit, ttb_min))
  }
  ttb_max <- ifelse(is.na(ttb_max),
    max(clinical_data_add$timeToTBFromVisit,
      na.rm = TRUE
    ),
    ttb_max
  )
  clinical_data_add <- clinical_data_add %>%
    dplyr::mutate(
      tfmttb = ifelse(Progressor == "no", # make zero non-progressor
        0,
        -9999
      ),
      tfmttb = ifelse(.data$Progressor == "yes", # make 0 if time to TB is greater than ttb_max; otherwise make
        # how many days it is from collapse time to the actual time
        pmax(0, .env$ttb_max - .data$timeToTBFromVisit),
        .data$tfmttb
      )
    )

  .data <- .data[, !(colnames(.data) %in%
    setdiff(
      colnames(clinical_data_add),
      cols_join
    ))]

  cn_vec_orig <- colnames(.data) %>%
    setdiff(cols_join)

  .data <- .data %>%
    dplyr::inner_join(
      clinical_data_add,
      by = cols_join
    )

  .data <- .data[, c(
    cols_join,
    cols_add,
    setdiff(
      cn_vec_orig,
      c(cols_join, cols_add)
    )
  ) %>%
    unique()]

  .data
}

#' @title Add assay data
add_tc_assay_data <- function(.data,
                              cols_add,
                              cols_join,
                              iter) {
  cn_orig <- colnames(.data)

  if (length(cols_add) == 0) {
    return(.data)
  }

  cols_add_split <- stringr::str_split(cols_add, "~")
  assay_vec <- purrr::map_chr(cols_add_split, function(x) x[1]) %>%
    unique()
  assay_to_cols_add <- purrr::map(assay_vec, function(x) {
    purrr::map_chr(cols_add_split, function(x) x[2])[
      purrr::map_chr(cols_add_split, function(x) x[1]) == x
    ]
  }) %>%
    setNames(assay_vec)

  for (assay in assay_vec) {
    .add_tc_assay_data <- switch(assay,
      "soma" = .add_tc_assay_data_soma,
      "risk6" = .add_tc_assay_data_risk6,
      "flow_ifng" = .add_tc_assay_data_flow_ifng,
      stop("assay not recognised")
    )
    .data <- .add_tc_assay_data(
      .data = .data,
      cols_add = assay_to_cols_add[[assay]],
      cols_join = cols_join,
      iter = iter
    )
  }

  .data[, c(
    cols_join,
    setdiff(colnames(.data), c(cn_orig)),
    setdiff(cn_orig, cols_join)
  )]
}

.add_tc_assay_data_flow_ifng <- function(.data,
                                         cols_add,
                                         cols_join,
                                         ...) {
  var_tbl_add <-
    TuberculomicsCompendium::tidy_wb_cd4p_infgp_cells_hladr_mfi

  var_tbl_add[, cols_add] <- var_tbl_add$CD4IFNg / var_tbl_add$CD4 * 1e2
  var_tbl_add <- var_tbl_add[, unique(c(cols_join, cols_add))]
  .data %>%
    dplyr::inner_join(
      var_tbl_add,
      by = cols_join
    )
}

.add_tc_assay_data_soma <- function(.data,
                                    cols_add,
                                    cols_join,
                                    ...) {
  var_tbl_add <- TuberculomicsCompendium::soma_data_tidy %>%
    dplyr::mutate(
      # prevents failures when fitting formula later
      Soma_Target = gsub("\\W", "", Soma_Target)
    ) %>%
    dplyr::filter(Soma_Target %in% cols_add) %>%
    dplyr::group_by_at(c(cols_join, "Soma_Target")) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  var_tbl_add <- var_tbl_add %>%
    tidyr::pivot_wider(
      id_cols = cols_join,
      names_from = "Soma_Target",
      values_from = "Soma_TransformedReadout"
    )

  var_tbl_add <- var_tbl_add[, c(cols_join, cols_add)]
  new_col_vec <- setdiff(colnames(var_tbl_add), cols_join)

  .data <- .data %>%
    dplyr::inner_join(
      var_tbl_add,
      by = cols_join
    )

  .data[, c(
    cols_join,
    new_col_vec,
    setdiff(
      colnames(.data),
      c(cols_join, new_col_vec)
    )
  )]
}

.add_tc_assay_data_risk6 <- function(.data,
                                     cols_add,
                                     cols_join,
                                     ...) {
  .data <- .data %>%
    dplyr::inner_join(
      TuberculomicsCompendium::signature_6gene,
      by = cols_join
    ) %>%
    dplyr::rename(risk6 = sig6gene_CorScore)
  .data[, c(cols_join, setdiff(colnames(.data), cols_join))]
}

.add_tc_assay_data_il2 <- function(.data,
                                   cols_add,
                                   cols_join,
                                   iter,
                                   ...) {
  data_join <- switch(cols_add,
    "il2_prob" = {
      DataTidyACSCyTOFCytokinesTCells::cd4_th1_il17$compass$locb0.15_min_clust[[iter$stim]]$fit
    },
    "il2_freq" = {

    },
    stop(paste0(paste0(cols_add, collapse = "/"), "not recognised"))
  )

  .data %>%
    dplyr::inner_join(
      data_join,
      by = cols_join
    )
}


#' @title Scale variables for modelling
scale_var <- function(.data, cols = NULL) {
  if (is.null(cols[1])) {
    return(.data)
  }

  if (length(cols) == 0) {
    return(.data)
  }

  if ("tfmttb" %in% cols && "tfmttb" %in% colnames(.data)) {
    .data <- .data %>%
      dplyr::mutate(
        origMeantfmttb = mean(tfmttb),
        tfmttb = tfmttb / 1e2
      )
  }
  if ("DaysSinceEntry" %in% cols && "DaysSinceEntry" %in% colnames(.data)) {
    .data <- .data %>%
      dplyr::mutate(DaysSinceEntry = (DaysSinceEntry - mean(DaysSinceEntry)) /
        sd(DaysSinceEntry))
  }
  .data
}


remove_tc_assay_from_exp_s <- function(p_dots) {
  nm_vec <- names(iter$var_exp_spline)
  for (i in seq_along(nm_vec)) {
    nm <- nm_vec[i]
    if (!grepl("^tc~", nm)) next
    nm_vec[i] <- stringr::str_split(nm, "~")[[1]][3]
  }
  names(iter$var_exp_spline) <- nm_vec
  p_dots
}

winsorise <- function(.data, wins, iter) {

  if (wins == "wins_n") {
    return(.data)
  }

  wins_var <- stringr::str_sub(wins, start = 6)
  wins_var_vec <- purrr::map_chr(
    seq_len(stringr::str_length(wins_var)), function(i) {
    stringr::str_sub(wins_var, i, i)
  })
  if ("y" %in% wins_var_vec) {
    var_vec_to_wins_y <- iter$var_dep
    for (i in seq_along(var_vec_to_wins_y)) {
      var <- var_vec_to_wins_y[[i]]
      if (!is.numeric(.data[[var]])) next
      num_vec <- .data[[var]]
      all_integer <- ceiling(num_vec) == num_vec & floor(num_vec) == num_vec
      if (!ind_list_null(iter$var_offset)) {
        num_vec_offsetted <- num_vec / .data[[iter$var_offset]]
        num_vec_wins <- .winsorise(x = num_vec_offsetted)
        num_vec_out <- num_vec_wins * .data[[iter$var_offset]]
      } else {
        num_vec_out <- .winsorise(x = num_vec)
      }
      if (all_integer) {
        num_vec_out <- round(num_vec_out)
      }
      .data[, var] <- num_vec_out
    }
  }
  if ("x" %in% wins_var_vec) {
    var_vec_to_wins_x <- c(iter$var_exp, .get_list_nm(iter$var_exp_spline))
    var_vec_to_wins_x <- .rm_list_null(var_vec_to_wins_x)
    for (i in seq_along(var_vec_to_wins_x)) {
      var <- var_vec_to_wins_x[[i]]
      if (!is.numeric(.data[[var]])) next
      num_vec <- .data[[var]]
      num_vec_out <- .winsorise(x = num_vec)
      .data[, var] <- num_vec_out
    }
  }
  .data
}



.winsorise <- function(x, mult = 3) {
  sd_var <- sd(x)
  mad_var <- mad(x)
  max_var <- max(
    mean(x) + mult * sd_var,
    median(x) + mult * mad_var
  )
  min_var <- min(
    mean(x) - mult * sd_var,
    median(x) - mult * mad_var
  )
  pmax(pmin(x, max_var), min_var)
}
