#' @title Install private GitHub package dependencies
#'
#' @param update logical.
#' If \code{FALSE}, then a package would
#' not be installed if it is already installed,
#' even if an update were available.
#' Default is \code{FALSE}.
#' @param project
#' "acs_cytof_tcells", "acs_cytof_nkbcells", "antibody".
#' Determines which packages are required for installation.
#' @param gitcreds_url
#' URL connected with GitHub PAT.
#' Passed to url argument of
#' \code{gitcreds::gitcreds_get}.
#' Default is \code{"https://github.com"}
install_pkg_private <- function(project,
                                update = FALSE,
                                gitcreds_url = "https://github.com" # nolint
                                ) {

  repo_vec_github_private <- switch(
    project,
    "acs_cytof_tcells" = c(
      "FredHutch/TuberculomicsCompendium",
      "SATVILab/pipeline",
      "SATVILab/PipelineAnalysisACS",
      "SATVILab/DataTidyACSCyTOFFAUST",
      "SATVILab/DataTidyACSCyTOFCytokinesTCells"
    ),
    "acs_cytof_nkbcells" = c(
      "FredHutch/TuberculomicsCompendium",
      "SATVILab/pipeline",
      "SATVILab/PipelineAnalysisACS",
      "SATVILab/DataTidyACSCyTOFFAUST",
      "SATVILab/DataTidyACSCyTOFCytokinesTCells"
    ),
    "antibody" = c(
      "FredHutch/TuberculomicsCompendium",
      "SATVILab/pipeline",
      "SATVILab/PipelineAnalysisACS",
      "SATVILab/DataTidyACSAntibody"
    ),
    stop(paste0("project ", project, " not recognised"), call. = FALSE)
  )

  repo_vec_github_private_to_install_ind <- vapply( # nolint
    repo_vec_github_private,
    function(x) {
      if (update) {
        return(TRUE)
      }
      if (is.null(names(x))) {
        pkg <- stringr::str_split(x, "/")[[1]][[2]]
      } else {
        pkg <- x
      }
      !(pkg %in% installed.packages())
    },
    logical(1)
  )
  repo_vec_github_private_to_install <- # nolint
    repo_vec_github_private[repo_vec_github_private_to_install_ind] # nolint
  if (length(repo_vec_github_private_to_install) > 0) {

    if (!requireNamespace("gitcreds", quietly = TRUE)) {
      install.packages("gitcreds")
    }
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github(
      repo_vec_github_private_to_install,
      auth_token = gitcreds::gitcreds_get(
        url = gitcreds_url,
        use_cache = TRUE
      )$password
    )
  }
  invisible(TRUE)
}
