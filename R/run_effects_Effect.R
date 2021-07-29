#' @title Run effects::effect when in an execution environment
#'
#' @return output from \code{effects::Effect(focal.predictors = var, mod = mod, xlevels = 50)}.
run_effects_Effect <- function(data, nm, var, mod) {

  if(exists(nm, envir = .GlobalEnv)){ # if there is already data_mod
    # in global env, then set data_mod back to that upon exit
    data_glbl_akc <- get(nm, envir = .GlobalEnv)
    on.exit(assign(nm, data_glbl_akc, envir = .GlobalEnv),
            add = TRUE)
  }

  assign(x = nm, value = data, envir = .GlobalEnv)


  # get estimates and ci's
  effects::Effect(
    focal.predictors = var,
    mod = mod,
    xlevels = 50)
}
