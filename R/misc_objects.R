#' @title Directory where data is processed
dir_data_prep <- "C:/Users/migue/Work/PhD/Code/cytofacs"

#' @title Get base directory to save results to
#'
#' @description Assumes that data are saved to file.path(pkg, "analysis", <sub_folder>), where the sub-folder's
#' name depends on the dataset and the purpose of the sub-folder (saving all results, saving "main" results for
#' sharing with colleagues or saving output directly for use in paper).
#'
#' Automatically creates folder if it does not yet exist.
#'
#' @param type "bulk", "share" or "paper". Specifies intention of saving. The values correspond to the folders
#' intended for saving all data, saving only data for main discussion with colleagues and saving data
#' for direction inclusion or processing for paper.
#' @param dataset character. Name of dataset being analysed.devt
#'
#' @details In future, could automatically add the resultant directory to .Rbuildignore and .gitignore (if not
#' already in it).
#' @param dataset "cytof". Specifies type of dataset.
get_dir_save <- function(type, dataset) {

  # package directory
  dir_proj <- DataPackageR::project_path()

  # directory to save to, based on type and dataset
  dir_save <- switch(type,
    "bulk" = file.path(dir_proj, "analysis", dataset),
    "share" = file.path(dir_proj, "analysis", paste0(dataset, "acs-analysis")),
    "paper" = file.path(dir_proj, "analysis", paste0(dataset, "acs-paper")),
    stop(paste0(type, " not one of bulk, share or paper in get_dir_save for base save directory"))
  )

  #
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  dir_save
}

#' @title Get directory where raw data are kept
#'
#' @description Gives directory where raw data, such as plots, are kept.
#' Assumes that the source folder is kept in same directory as acsanalysis package,
#' and that the final folder in the directory is paste0(dataset, "acs").
#'
#' @param dataset "cytof" or "antibody".
#'
#' @return Character.
get_dir_data_prep <- function(dataset) {
  if (!dataset %in% c("cytof", "antibody")) {
    stop(paste0(dataset, " (dataset parameter in get_dir_data_prep fn) should be one of cytof or antibody"))
  }
  dir_analysis <- DataPackageR::project_path()
  slash_loc_tbl <- stringr::str_locate_all(dir_analysis, "/")[[1]]
  last_slash_loc <- slash_loc_tbl[nrow(slash_loc_tbl), "end"][[1]]
  dir_base <- str_sub(dir_analysis, end = last_slash_loc - 1)
  dir_data_prep <- file.path(dir_base, paste0(dataset, "acs"), "data-raw")
  if (dataset == "cytof") dir_data_prep <- file.path(dir_data_prep, "gating")
  dir_data_prep
}
