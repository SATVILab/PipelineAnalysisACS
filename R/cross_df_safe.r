#' @title Safely cross lists
#'
#' @description Returns the cartesian of all its elements in a list
#' as a dataframe, with one combination per row.

#'
#' @param .l A list of lists or atomic vectors (or, a dataframe).
#' @param .filter function.
#' A predicate function that takes the same number of arguments as the number of variables to be combined.
#'
#' @details
#' 
#' Vendored from `MiguelRodo/datautils`.
#' 
#' Elements in the list to be crossed
#' such as \code{list(1, 2:3)} keep
#' their multi-length elements in one combination
#' (i.e. one dataframe row).
#' That is the sense in which is "safe",
#' as it behaves as expected in that
#' not only lists are "guarded", i.e. not joined together.
#' 
#' The function `purrr::cross` returns an error
#' in this case, whilst `purrr::cross` piped
#' `purrr::map_df(tibble::as_tibble)` w
#' creates two combinations (two dataframe rows)
#' from the second list element above.
#' One could wrap the second element in a list, but
#' that adds an extra layer when accessing the element
#' in the combination.
#'
#' Another advantage
#'
#' See examples for a demonstration.
#'
#' @examples
#' 
#' # cross_df_safe automatically
#' # "protects" non-list list elements
#' to_cross_list <- list(
#'  V1 = c("a", "b"),
#'  V2 = list(3, 4:5)
#' )
#' # returns an error
#' try(purrr::cross_df(to_cross_list))
#' # can avoid this by wrapping each to-cross
#' # element's elements in lists,
#' # but then the crossed dataframe
#' # has second column's elements
#' # each be lists
#' to_cross_list_list <- list(
#'   V1 = c("a", "b"),
#'   V2 = list(list(3), list(4:5))
#' )
#' cross_df_purrr <- purrr::cross_df(to_cross_list_list)
#' cross_df_purrr$V2
#' # returns a four-row dataframe, where
#' # second column has elements of lengths 1, 1, 2 and 2
#' # that are not lists
#' cross_df_sf <- cross_df_safe(to_cross_list)
#' cross_df_sf$V2
#' # the above is more natural to work with
#' # no element in a combination
#' # is needlessly wrapped inside an unnamed list of length 1
#' 
#' # 
#' to_cross_list_list <- list(
#'   V1 = c(list("a"), list("b")),
#'   V2 = list(
#'     c("a" = "c"),
#'     c("b" = "d")
#'  )
#' )
#'
#'cross_df_std <- purrr::cross_df(to_cross_list_list)
#'cross_df_sf <- cross_df_safe(to_cross_list_list)
#' # purr:cross_df needlessly places each element of V2
#' # (e.g. list("a" = "c")) in its own list:
#' cross_df_std$V2
#' #the above means that to refer to the element, you need to 
#' # add an extra set of brackets:
#' # cross_df_std$V2[[1]][[1]]
#' # in this case, however,cross_df_safe flattens the list by one level:
#' cross_df_sf$V2
#' # so that you can refer to individual elements more succinctly:
#' cross_df_sf$V2[[1]]
#' @return A dataframe with the same number of rows as
#'
#' @export
cross_df_safe <- function(.l, .filter = NULL) {
  cols_to_protect_vec <- get_cols_to_protect(.l)

  cross_list <- purrr::cross(.l, .filter)

  cross_df <- cross_list |>
    purrr::map_df(.convert_to_tibble_row_one,
                  protect = cols_to_protect_vec)

  list_col_ind <- which(purrr::map_lgl(cross_df, is.list))

  for (i in seq_along(list_col_ind)) {
    ind <- list_col_ind[[i]]
    all_list_of_list <- purrr::map_lgl(
      cross_df[[ind]],
      is.list
    ) |>
      all()
    all_unnamed <- is.null(names(cross_df[[ind]]))
    all_length_1 <- purrr::map_lgl(
      cross_df[[ind]],
      function(x) length(x) == 1
    ) |>
      all()

    if (all_list_of_list && all_unnamed && all_length_1) {
      cross_df[, ind] <- list(purrr::flatten(cross_df[[ind]]))
    }
  }

  cross_df
}

get_cols_to_protect <- function(.l) {
  purrr::map_lgl(
    seq_along(.l),
    function(ind) {
      list_vec_lgl <- purrr::map_lgl(.l[[ind]], function(x) {
        is.list(x)
      })
      if (any(list_vec_lgl)) {
        return(TRUE)
      }
      elem_vec_lgl <- purrr::map_lgl(.l[[ind]], function(x) {
        is.numeric(x) || is.logical(x) || is.character(x)
      })
      if(!all(elem_vec_lgl)) return(FALSE)

      lg1 <- purrr::map_lgl(.l[[ind]], function(x) {
        length(x) > 1
      }) |>
        any()

      if (any(lg1)) {
        return(TRUE)
      }
      named <- purrr::map_lgl(.l[[ind]], function(x) {
        !is.null(names(x))
      }) |>
        all()
      named
    }
  ) |>
    which()
}

.convert_to_tibble_row_one <- function(x, protect = NULL) {
  for (i in seq_along(protect)) {
    x[[protect[i]]] <- list(x[[protect[i]]])
  }
  tibble::as_tibble(x)
}
