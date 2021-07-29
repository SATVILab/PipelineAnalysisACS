#' @title Make automatic plot using effects package
#'
#' @param mod model object.
#' @param data_mod dataframe.
#' @param var character vector. Names of variables to be saved. Only required
#' if height or width is NULL, as its length is used to calculate
#' them automatically.
#' @param fn character. File name. Default is \code{"p_fit_auto.png"}.
#' @param dir_save character. Folder to save to.
#' @param height,width numeric. If not \code{NULL}, then
#' the plot height and width are given by these.
#' Default is \code{NULL}.
#'
#' @return Returns path to save to invisibly.
#'
#' @details
#' Note that only PNG output is supported at present.
#' The plot function for `effects::predictorEffects` assumes
#' that the data used to fit the model
#' are in the global environment. This function
#' automatically assigns the data to the global environment
#' to an object with the correct name, and restores
#' any object "overwritten" after the function is completed running.
#'
#' @export
plot_display_auto <- function(mod, data_mod, var, fn, dir_save, height = NULL, width = NULL) {

  if (is.null(height) || is.null(width)) {
    n_var <- length(var)
  }
  if (is.null(height)) height <- ceiling(n_var / 2) * 8
  if (is.null(width)) width <- ifelse(n_var >= 2, 19, 9.5)

  png(filename = file.path(dir_save, "p_fit_auto.png"), units = 'cm',
      height = height, width = width, res = 300,
      pointsize = 10)

  # add data_mod to global env and ensure removal upon function exit
  # --------------------------
  on.exit(suppressWarnings(try(
    rm(data_mod, envir = .GlobalEnv),
    silent = TRUE)))

  if(exists('data_mod', envir = .GlobalEnv)){ # if there is already data_mod
    # in global env, then set data_mod back to that upon exit
    data_mod_old_akx <- get('data_mod', envir = .GlobalEnv)
    on.exit(assign('data_mod', data_mod_old_akx, envir = .GlobalEnv),
            add = TRUE)
  }
  assign(x = "data_mod", value = data_mod, envir = .GlobalEnv)

  # create standard predictorEffects plot
  plot_out <- try(plot(effects::predictorEffects(mod = mod)))
  invisible(dev.off())
  if(identical(class(plot_out), 'try-error')){
    file.remove(file.path(dir_save, "p_fit_auto.png"))
  }
  rm(plot_out)
  invisible(file.path(dir_save, "p_fit_auto.png"))
}
