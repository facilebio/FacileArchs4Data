# An S3 wrapper to Arsch4data
# (trying to use R6 makes integration into rest of ecosystem more difficult)

#' Constructor
#' @export
#' @param x name of species, of path to the h5 file
#' @examples
#' mfds <- FacileArchs4DataSet("mouse")
#'
FacileArchs4DataSet <- function(
  x,
  ...,
  species = NULL,
  data_dir = getOption("archs4.data_dir", NULL),
  use_cache = TRUE,
  remove_sc = TRUE,
  threshold_sc = 0.5
) {
  if (!test_file_exists(x, extension = c("h5", "hdf5"))) {
    x <- match.arg(x, c("human", "mouse"))
    species <- x
    assert_directory_exists(data_dir, "r")
    fn.regex <- sprintf(".*%s.*\\.(h5|hdf5)$", species)
    x <- dir(data_dir, fn.regex, full.names = TRUE)
    if (length(x) == 0L) {
      stop("Can not archs4 hdf5 file for species: ", species)
    }
    if (length(x) > 1L) {
      stop(
        "More than one hdf5 file found for speies: ",
        paste(x, collapse = ",")
      )
    }
  }
  assert_file_exists(x, "r")
  if (is.null(species)) {
    if (grepl("mouse", basename(x), ignore.case = TRUE)) {
      species <- "mouse"
    } else {
      species <- "human"
    }
  }

  samples.all <- .load_archs4_samples(x, use_cache = TRUE) |>
    dplyr::rename(dataset = "series_id", sample_id = "sample")
  class(samples.all) <- c(
    "archs4_facile_frame",
    "facile_frame",
    class(samples.all)
  )

  out <- list(
    h5 = x,
    samples = samples.all,
    features = .load_archs4_features(x, species),
    species = species,
    remove_sc = assert_flag(remove_sc),
    threshold_sc = assert_number(threshold_sc, lower = 0, upper = 1)
  )

  out$sample_source <- attr(out$samples, "source")

  class(out) <- c("FacileArchs4DataSet", "FacileDataStore", class(out))
  out$study_distro <- dplyr::count(archs4_studies(out), likely_sc)
  out$samples <- FacileData::set_fds(out$samples, out)
  out
}

#' Enumerate the studies in an archs4 dataset
#'
#' @export
#' @param x the FacileArchs4DataSet
#' @param sample_search If a search string is provided, this will trigger
#'   a sample search via the [archs4_sample_search()] method. The study table
#'   returned will be augmented with a `search_hits` column that tells you how
#'   many samples within the study were a hit for your search.
#' @param scrub_search,meta_fields,ignore.case parameters passed to `search` function.
#' @param threshold_sc The singlecellprobability threshold to break for a sample
#'   to be considered likely single cell. This defaults to `x$threshold_sc`
#' @param threshold_sc_likely the proportion of samples in a study that
#'   must look like singlecell data in order to classify the study as
#'   singlecell. default: 0.2
archs4_studies <- function(
  x,
  sample_search = NULL,
  meta_fields = .default_meta_fields(),
  threshold_sc = x$threshold_sc,
  threshold_sc_likely = 0.2,
  ...
) {
  assert_class(x, "FacileArchs4DataSet")
  assert_number(threshold_sc, lower = 0.01, upper = 1)
  assert_number(threshold_sc_likely, lower = 0.05, upper = 1)

  hits <- NULL
  if (checkmate::test_string(sample_search)) {
    hits <- archs4_sample_search(x, sample_search, meta_fields, ignore.case) |>
      dplyr::summarize(search_hits = dplyr::n(), .by = "series_id")
  }
  smry <- samples(x) |>
    dplyr::summarize(
      nsamples = dplyr::n(),
      nsc = sum(singlecellprobability >= threshold_sc),
      likely_sc = nsc / nsamples >= threshold_sc_likely,
      search_hits = 0L,
      .by = "dataset"
    )
  if (!is.null(hits) && nrow(hits) > 0L) {
    smry <- smry |>
      dplyr::select(-search_hits) |>
      dplyr::left_join(hits, by = "dataset") |>
      dplyr::mutate(search_hits = ifelse(is.na(search_hits), 0L, search_hits))
  }
  smry
}

#' Find samples that match a regular expression in any of the descriptive
#' columns.
#'
#' This function is fashioned after the archs4py.meta() function.
#'
#' @seealso https://github.com/MaayanLab/archs4py/blob/main/archs4py/data.py#L27
#' @param x A FacileArchs4DataSet or a subset of its sample frame
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
#' @return the samples facile_frame that matches the search
archs4_sample_search <- function(
  x,
  regex,
  meta_fields = .default_meta_fields(),
  scrub_search = FALSE,
  ignore.case = TRUE,
  remove_sc = x$remove_sc,
  threshold_sc = x$threshold_sc
) {
  xdf <- if (is(x, "FacileArchs4DataSet")) samples(x) else x
  assert_data_frame(xdf)

  out <- fuzzy_search(
    xdf,
    regex = regex,
    meta_fields = meta_fields,
    scrub_search = scrub_search,
    ignore.case = ignore.case,
    remove_sc = remove_sc,
    threshold_sc = threshold_sc
  )
  out
}

#' @noRd
#' @export
print.FacileArchs4DataSet <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
format.FacileArchs4DataSet <- function(x, ...) {
  ntotal <- sum(x$study_distro$n)
  nbulk <- sum(x$study_distro$n[!x$study_distro$likely_sc])
  nsc <- ntotal - nbulk
  out <- paste(
    "=======================================================================\n",
    name(x),
    "\n",
    "-----------------------------------------------------------------------\n",
    sprintf("  path: %s\n", x$h5),
    sprintf("  studies: %s\n", prettyNum(ntotal, big.mark = ",")),
    sprintf("  ├─ bulk: %s\n", prettyNum(nbulk, big.mark = ",")),
    sprintf("  └─ single cell: %s\n", prettyNum(nsc, big.mark = ",")),
    sprintf("  samples: %s\n", prettyNum(nrow(x$samples), big.mark = ",")),
    sprintf("  genes: %s\n", prettyNum(nrow(x$features), big.mark = ",")),
    "  .....................................................................\n",
    sprintf("  sample source: %s\n", x$sample_source),
    "=======================================================================\n",
    sep = ""
  )
}
