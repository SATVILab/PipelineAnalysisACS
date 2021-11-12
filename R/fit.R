#' @import lme4
#' @import glmmTMB
#'
#' @export
fit <- function(data_mod, data_raw, p_dots, dir_proj, iter) {
  if (nrow(data_mod) == 0) {
    stop("0 rows in data_mod")
  }

  # ============================
  # Preparation
  # ============================

  # list to save results to
  # ---------------------------
  mod_list <- list()

  # Model formulae
  # ---------------------------

  if (identical(iter[["var_offset"]], "none") ||
    !"var_offset" %in% colnames(iter)) {
    iter$var_offset <- NULL
  }
  if (identical(iter[["var_conf"]], "none") ||
    !"var_conf" %in% colnames(iter)) {
    iter$var_conf <- NULL
  }
  fml_as_chr_list <- modutils::get_mod_fml_as_chr(
    var_dep = "resp",
    var_exp = iter[["var_exp"]] %>% unlist(),
    var_exp_spline = iter[["var_exp_spline"]],
    var_re = iter[["var_re"]] %>% unlist(),
    var_offset = switch(as.character(is.null(iter[["var_offset"]])),
      "TRUE" = NULL,
      "FALSE" = paste0("log(", iter$var_offset, ")")
    ),
    var_conf = iter[["var_conf"]] %>% unlist(),
    var_int = iter[["var_int"]] %>% unlist(),
    rhs_text = NULL,
    int = TRUE,
    sub_exp = iter$sub_exp
  )

  # ============================
  # Fit models
  # ============================

  mod_fn <- switch(iter$pkg,
    "lme4" = switch(iter$family,
      "nb" = paste0("lme4::glmer.nb"),
      "normal" = paste0("lme4::lmer"),
      "lme4::glmer"
    ),
    "glmmTMB" = "glmmTMB::glmmTMB",
    paste0(iter$pkg, "::", iter$family)
  )


  family_arg <- switch(as.character(is.null(iter$family)),
    "FALSE" = switch(iter$pkg,
      "lme4" = switch(iter$family,
        "normal" = ,
        "nb" = NULL,
        "Gamma" = ,
        "gamma" = "Gamma(link = 'log')",
        paste0("family = ", iter$family)
      ),
      "glmmTMB" = paste0(
        "family = ",
        switch(iter$family,
          "normal" = "gaussian()",
          "Beta" = ,
          "beta" = "beta_family()",
          "gamma" = ,
          "betabin" = "betabinomial(link = 'logit')",
          "bin" = "binomial",
          "nb" = "nbinom2",
          "Gamma" = "Gamma(link = 'log')",
          stop("family does not match one of the available ones for glmmTMB")
        )
      ),
      stop("pkg not one of the available pkg's")
    ),
    "TRUE" = NULL
  )

  family_arg <- switch(as.character(iter$family %in% c("betabin", "bin")),
    "TRUE" = paste0(family_arg, ", weights = data_mod$n_cell"),
    "FALSE" = family_arg
  )

  mod_arg_txt <- list(
    family_arg,
    iter[["mod_args"]],
    "data = data_mod"
  ) %>%
    purrr::compact() %>%
    unlist() %>%
    paste0(collapse = ", ")

  # fit full model
  # ----------------

  # fit it initially
  mod_expr_as_chr_full <- paste0(
    "try(",
    mod_fn, "(",
    fml_as_chr_list$full, ", ",
    mod_arg_txt,
    "))"
  )

  mod_full <- eval(rlang::parse_expr(mod_expr_as_chr_full))


  if (class(mod_full) == "try-error") {
    mod_full <- modutils::refit_mod(
      mod_expr_as_chr = mod_expr_as_chr_full,
      pkg = iter$pkg,
      family = iter$family
    )
  }

  if (class(mod_full) == "try-error") {
    return(list(full = mod_full))
  }

  # get theta of full model for use in later models
  theta_fix <- switch(as.character(iter$family == "nb" && iter$pkg == "lme4"),
    "TRUE" = lme4::getME(mod_full, "glmer.nb.theta"),
    "FALSE" = NULL,
    stop("switch statement incorrectly specified")
  )

  # fit remaining models
  non_full_fml_ind_vec <- which(names(fml_as_chr_list) != "full")
  mod_list_non_full <- purrr::map(fml_as_chr_list[non_full_fml_ind_vec], function(fml_as_chr) {
    mod_expr_as_chr <- paste0(
      "try(",
      mod_fn, "(",
      fml_as_chr, ", ",
      mod_arg_txt,
      "))"
    )
    mod <- eval(rlang::parse_expr(mod_expr_as_chr))
    if (class(mod) == "try-error") {
      mod_fix <- modutils::refit_mod(
        mod_expr_as_chr = mod_expr_as_chr,
        pkg = iter$pkg,
        family = iter$family,
        theta_fix = theta_fix
      )
    } else {
      return(mod)
    }
    if (class(mod_fix) == "try-error" && class(mod) == "try-error") {
      return(NULL)
    }
    mod_fix
  }) %>%
    setNames(names(fml_as_chr_list)[non_full_fml_ind_vec])

  mod_list <- list("full" = mod_full) %>%
    append(mod_list_non_full)

  # ============================
  # Save
  # ============================

  save_objects(
    obj_list = list("mod_list" = mod_list),
    dir_proj = iter$dir_stg
  )


  mod_list
}
