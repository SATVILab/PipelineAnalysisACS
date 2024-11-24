% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc_objects.R
\name{get_dir_save}
\alias{get_dir_save}
\title{Get base directory to save results to}
\usage{
get_dir_save(type, dataset)
}
\arguments{
\item{type}{"bulk", "share" or "paper". Specifies intention of saving. The values correspond to the folders
intended for saving all data, saving only data for main discussion with colleagues and saving data
for direction inclusion or processing for paper.}

\item{dataset}{"cytof". Specifies type of dataset.}
}
\description{
Assumes that data are saved to file.path(pkg, "analysis", <sub_folder>), where the sub-folder's
name depends on the dataset and the purpose of the sub-folder (saving all results, saving "main" results for
sharing with colleagues or saving output directly for use in paper).

Automatically creates folder if it does not yet exist.
}
\details{
In future, could automatically add the resultant directory to .Rbuildignore and .gitignore (if not
already in it).
}
