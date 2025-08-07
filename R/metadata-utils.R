#' The `characteristics_ch1` column often has key:value attributes worth parsing.
#'
#' Often times the archs4 `characteristics_ch1` sample column looks like
#' `"cell line: K562,cell type: Immortalized epithelial cells,genotype: WT, ...`
#' Let's bust this up into usable data.frame columns.
#'
#' @export
#' @param x the samples archs4 data.frame
#' @param column the name of the column to key:value extract, default
#'   is `"characteristics_ch1"`
#' @param delim_var the delimeter used to split key:value pairs, default: `","`
#' @param delim_key the delimeter to split key/value pairs, default: `":"`
#' @param clean_names use [janitor::clean_names()] to clean up the key variable
#'   names expanded from `column`.
#' @param keep_column return the `x[[column]]` column in the output, defaults
#'   to `TRUE`
#' @param ... stuff and things
#' @return wide tibble of metadata variables from `x[[column]]`
#' @examples
#' a4 <- FacileArchs4Data::FacileArchs4DataSet("human")
#' info <- samples(a4) |> dplyr::sample_n(20) |> arcsh4_split_metadata()
arcsh4_split_metadata <- function(
  x,
  column = "characteristics_ch1",
  delim_var = ",",
  delim_key = ":",
  clean_names = TRUE,
  keep_column = TRUE,
  ...
) {
  if (FALSE) {
    a4 <- FacileArchs4Data::FacileArchs4DataSet("human")
    x <- samples(a4) |> dplyr::sample_n(20)
    column <- "characteristics_ch1"
    delim_var <- ","
    delim_key <- ":"
    keep_column <- TRUE
  }
  checkmate::assert_character(x[["dataset"]])
  checkmate::assert_character(x[["sample_id"]])
  checkmate::assert_character(x[[column]])
  info <- x |>
    dplyr::select(dataset, sample_id, dplyr::all_of(column)) |>
    tidyr::separate_longer_delim(column, delim = delim_var) |>
    tidyr::separate_wider_delim(
      column,
      delim = delim_key,
      names = c("key", "value"),
      too_many = "merge",
      cols_remove =!keep_column
    ) |>
    dplyr::mutate(
      key = janitor::make_clean_names(key, allow_dupes = TRUE),
      value = trimws(value)) |>
    tidyr::pivot_wider(names_from = "key", values_from = "value")
  if (keep_column) {
    info <- info |>
      dplyr::relocate(dplyr::all_of(column), .after = dplyr::last_col())
  }
}
