#' @title Open a dialogue box at selected folder
#'
#' @description Useful for easily opening a window in
#' which you can see all the files in a folder and open them using
#' \code{Open with...}. The default is to open the folder
#' at the specific project directory for the present combination
#' of tuning parameters.
#'
#' @param path character. The initial working directory, from which the file dialog should begin browsing. Defaults to the current RStudio project directory.
#'
#' @inheritDotParams rstudioapi::selectFile
#'
#' @return \code{invisible(<selected_file>/NULL)} (\code{NULL} if no file were
#'   selected).
#'
#' @details
#' Only works in RStudio. Wrapper around `rstudioapi::selectFile` that simply
#' passes all its arguments to it but avoids printing the path to a file
#' if a file were selcted.
#'
#' @export
#' @aliases open_dir
#' @examples
#' open_folder()
#' open_dir()
#' open_folder(.libPaths()[1])
open_folder <- function(path = parent.frame(2)$dir_proj, ...) invisible(rstudioapi::selectFile(path = path, ...))

#' @rdname open_folder
open_dir <- function(path = parent.frame(2)$dir_proj,...) invisible(rstudioapi::selectFile(path = path, ...))

#' @rdname open_folder
open_proj <- function(path = parent.frame(2)$dir_proj,...) invisible(rstudioapi::selectFile(path = path, ...))
