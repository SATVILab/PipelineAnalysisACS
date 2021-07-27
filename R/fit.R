#' @import lme4
#' @import glmmTMB
#'
#' @export
fit <- function(data_mod, data_raw, p_dots, dir_proj){

  # ============================
  # Preparation
  # ============================

  # list to save results to
  # ---------------------------
  mod_list <- list()

  # Model formulae
  # ---------------------------

  if(identical(p_dots$var_offset, 'none')) p_dots$var_offset <- NULL
  if(identical(p_dots$var_conf, 'none')) p_dots$var_conf <- NULL

  fml_as_chr_list <- modutils::get_mod_fml_as_chr(var_dep = p_dots$var_dep,
                                                  var_exp = p_dots$var_exp,
                                                  var_exp_spline = p_dots$var_exp_spline,
                                                  var_re = p_dots$var_re,
                                                  var_offset = switch(as.character(is.null(p_dots$var_offset)),
                                                                      "TRUE" = NULL,
                                                                      "FALSE" = paste0("log(", p_dots$var_offset, ")")),
                                                  var_conf = p_dots$var_conf,
                                                  var_int = p_dots$var_int,
                                                  rhs_text = NULL,
                                                  int = TRUE,
                                                  sub_exp = p_dots$sub_exp)

  # ============================
  # Fit models
  # ============================

  mod_fn <- switch(p_dots$pkg,
                   "lme4" = switch(p_dots$family,
                                   "nb" = paste0("lme4::glmer.nb"),
                                   "normal" = paste0("lme4::lmer"),
                                   "lme4::glmer"),
                   "glmmTMB" = "glmmTMB::glmmTMB",
                   paste0(p_dots$pkg, "::", p_dots$family))


  family_arg <- switch(as.character(is.null(p_dots$family)),
                       "FALSE" = switch(p_dots$pkg,
                                        'lme4' = switch(p_dots$family,
                                                        "normal" = ,
                                                        "nb" = NULL,
                                                        "Gamma" = ,
                                                        "gamma" = "Gamma(link = 'log')",
                                                        paste0('family = ', p_dots$family)),
                                        'glmmTMB' = paste0('family = ',
                                                           switch(p_dots$family,
                                                                  'beta' = 'beta_family()',
                                                                  'gamma' = ,
                                                                  'Gamma' = "Gamma(link = 'log')",
                                                                  stop("family does not match one of the available ones for glmmTMB"))),
                                        stop("pkg not one of the available pkg's")),
                       "TRUE" = NULL)

  mod_arg_txt <- list(
    family_arg,
    p_dots$mod_args,
    'data = data_mod') %>%
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


  if(class(mod_full) == 'try-error'){
    mod_full <- modutils::refit_mod(mod_expr_as_chr = mod_expr_as_chr_full,
                                    pkg = p_dots$pkg,
                                    family = p_dots$family)
  }

  if(class(mod_full) == 'try-error') return(list(full = mod_full))

  # get theta of full model for use in later models
  theta_fix <- switch(as.character(p_dots$family == 'nb' && p_dots$pkg == 'lme4'),
                      "TRUE" = lme4::getME(mod_full, "glmer.nb.theta"),
                      "FALSE" = NULL,
                      stop("switch statement incorrectly specified"))

  # fit remaining models
  non_full_fml_ind_vec <- which(names(fml_as_chr_list) != 'full')
  mod_list_non_full <- purrr::map(fml_as_chr_list[non_full_fml_ind_vec], function(fml_as_chr){
    mod_expr_as_chr <- paste0(
      "try(",
      mod_fn, "(",
      fml_as_chr, ", ",
      mod_arg_txt,
      "))"
    )
    mod <- eval(rlang::parse_expr(mod_expr_as_chr))
    if(class(mod) == 'try-error'){
      mod_fix <- modutils::refit_mod(mod_expr_as_chr = mod_expr_as_chr,
                                     pkg = p_dots$pkg,
                                     family = p_dots$family,
                                     theta_fix = theta_fix)
    } else return(mod)
    if(class(mod_fix) == 'try-error' && class(mod) == 'try-error'){
      return(NULL)
    }
    mod_fix
  }) %>%
    setNames(names(fml_as_chr_list)[non_full_fml_ind_vec])

  mod_list <- list('full' = mod_full) %>%
    append(mod_list_non_full)

  # ============================
  # Save
  # ============================

  save_objects(obj_list = list('mod_list' = mod_list),
               dir_proj = p_dots$dir_stg)


  mod_list
}

