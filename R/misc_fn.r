.ind_list_null <- function(x) {
  vapply(x, function(.x) identical(list(NULL), .x), logical(1))
}

.rm_list_null <- function(x) {
  x[!.ind_list_null(x)]
}

.get_list_nm <- function(x) {
  lapply(x, function(.x) names(.x))
}