#' Fuzzy search across columns to find samples that match search term
#'
#' Search for samples that are "kind of like" soething you are interested in
#' @export
#' @param x a tibble of sample information
#' @param regex the regular expression to use for search
#' @param meta_fields which columns in the dataset to search against, defaults
#'   to a set of columns defined in the archs4py package
#' @param remove_sc Wether to remove single cell samples from search?
#'   Default: `TRUE`
#' @param ... stuff and things
#' @param scrub_search remove whitespace,dashes,etc. from values in
#'   samples columns before running the regex? Setting this to `TRUE` slows
#'   down the search appreciable (~ 3x - 4x). Default: `FALSE`
#' @param ignore.case Should we ignore case when using `regex`?
#'   Default: `TRUE`
#' @param threshold_sc the probability threshold to indicate wheter sample
#'   is from a singlecell dataset, default: 0.5
#' @return the samples tibble that matches the search
fuzzy_search <- function(
  x,
  regex,
  meta_fields = NULL,
  scrub_search = FALSE,
  ignore.case = TRUE,
  remove_sc = TRUE,
  threshold_sc = 0.5
) {
  assert_data_frame(x, "data.frame")
  assert_string(regex)
  assert_flag(remove_sc)
  if (is.null(meta_fields)) {
    meta_fields <- .default_meta_fields()
  }
  assert_character(meta_fields)
  mfields <- intersect(meta_fields, colnames(x))
  assert_character(mfields, min.len = 1L)
  bad.fields <- setdiff(meta_fields, mfields)
  if (length(bad.fields) > 0) {
    warning(
      "These search fields do not exist and will be ignored: ",
      paste(bad.fields, collapse = ",")
    )
  }
  assert_flag(scrub_search)
  assert_flag(remove_sc)
  assert_flag(ignore.case)
  assert_number(threshold_sc, lower = 0.01, upper = 1)
  if (remove_sc) {
    samples <- dplyr::filter(x, .data$singlecellprobability < threshold_sc)
  } else {
    samples <- x
  }
  hits <- rep(FALSE, nrow(samples))
  for (cname in meta_fields) {
    text <- samples[[cname]]
    if (scrub_search) {
      text <- gsub("_|-|'|/| |\\.", "", text)
    }
    hits <- hits | grepl(regex, text, ignore.case = TRUE)
  }
  samples[hits, , drop = FALSE]
}
