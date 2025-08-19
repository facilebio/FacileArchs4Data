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
  threshold_sc = 0.5,
  .duckdb = FALSE
) {
  if (!test_file_exists(x, extension = c("h5", "hdf5"))) {
    x <- match.arg(x, c("human", "mouse"))
    species <- x
    info <- archs4_dir_info(data_dir) |>
      dplyr::filter(.data$species == .env$species)
    if (nrow(info) == 0L) {
      stop("Can not archs4 hdf5 file for species: ", species)
    }
    if (nrow(info) > 1L) {
      warning("Multiple files found for species, taking latest version")
      info <- info[1,]
    }
    out <- list(
      h5 = info$path,
      path.sample_cache = unname(info$sample_cache_path),
      species = info$species
    )
  } else {
    out <- list(
      h5 = x,
      path.sample_cache = unname(archs4_cache_fn(x, outdir = dirname(x))),
      species = .infer_species_from_filename(x)
    )
  }

  out$path.series_info <- file.path(
    sub("samples.parquet$", "cache", out$path.sample_cache),
    "series_info"
  )

  assert_file_exists(out$h5, "r")
  out$meta <- .load_archs4_metadata(out$h5)
  out$remove_sc <- assert_flag(remove_sc)
  out$threshold_sc = assert_number(threshold_sc, lower = 0, upper = 1)
  out$features <- .load_archs4_features(out$h5, out$species)

  samples.all <- .load_archs4_samples(
    out$h5,
    use_cache = TRUE,
    sample_cache_path = out$path.sample_cache,
    create_cache = create_cache,
    .duckdb = .duckdb,
    ...
  )

  out$samples <- samples.all
  out$duckdb <- !test_data_frame(samples.all)
  out$sample_source <- attr(out$samples, "source")
  class(out) <- c("FacileArchs4DataSet", "FacileDataStore", class(out))
  out$study_distro <- dplyr::count(archs4_studies(out), likely_sc)

  class(out$samples) <- c(
    "archs4_facile_frame",
    "facile_frame",
    class(out$samples)
  )
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
  if (checkmate::test_string(sample_search) && !x$duckdb) {
    hits <- archs4_sample_search(x, sample_search, meta_fields, ignore.case) |>
      dplyr::summarize(search_hits = dplyr::n(), .by = "series_id")
  }

  if (x$duckdb) {
    smry <- samples(x) |>
      dplyr::summarize(
        nsamples = dplyr::n(),
        nsc = sum(singlecellprobability >= threshold_sc, na.rm = TRUE),
        .by = "dataset"
      ) |>
      dplyr::collect() |>
      dplyr::mutate(
        likely_sc = nsc / nsamples >= threshold_sc_likely,
        search_hits = 0L,
      )
  } else {
    smry <- samples(x) |>
      dplyr::summarize(
        nsamples = dplyr::n(),
        nsc = sum(singlecellprobability >= threshold_sc, na.rm = TRUE),
        likely_sc = nsc / nsamples >= threshold_sc_likely,
        search_hits = 0L,
        .by = "dataset"
      )
  }

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
#' @export
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

#' Extract metadata from arcsh4 dataset
#'
#' @export
#' @param x a FacileArchs4DataSet or the path to the source hdf5 data file
#' @param variable the name of the variable to extract from the metadata field.
#'   If `NULL` (default), all values will be returned.
#' @param ... stuff and things
#' @param multiple_ok If `TRUE`, then `x` can be a pointer to many hdf5 files
#' @return a list of metadata pulled from the archs4 hdf5 file. If `x` is
#'   a list of multiple hdf5 files, then a tibble of metadata is returned,
#'   each row correxponds to the file passed in `x`
archs4_meta <- function(
  x,
  variable = NULL,
  ...,
  multiple_ok = is.character(x)
) {
  assert_character(variable, null.ok = TRUE)
  if (is(x, "FacileArchs4DataSet")) {
    meta <- list(x$meta)
  } else {
    assert_character(x, min.len = 1L, max.len = if (multiple_ok) NULL else 1L)
    assert_file_exists(x, extension = "h5")
    meta <- lapply(x, .load_archs4_metadata)
  }
  if (is.null(variable) || length(variable) == 0L) {
    variable <- names(meta[[1L]])
  }
  var.bad <- setdiff(variable, names(meta[[1L]]))
  if (length(var.bad)) {
    warning("Unknown meta variables: ", paste(var.bad, collapse = ","))
    variable <- setdiff(variable, var.bad)
  }
  meta <- lapply(meta, "[", variable)
  if (length(meta) == 1L) meta[[1L]] else dplyr::bind_rows(meta)
}

#' Return the version of an arcsh4 dataset
#' @export
#' @inheritParams archs4_meta
#' @return the version string from the archs4 hdf5 file
archs4_version <- function(x, ...) {
  archs4_meta(x)$version
}

#' Return the registered creation date of an arcsh4 dataset
#' @export
#' @inheritParams archs4_meta
#' @param string return the data as a string (default) or a `Date` object
#' @param ... stuff and things
#' @param format the format return the date in
#' @return the creation date registered in archs4 hdf5 file
archs4_creation_date <- function(x, string = TRUE, ..., format = "%Y-%m-%d") {
  assert_flag(string)
  assert_string(format)
  out <- archs4_meta(x)$creation_date
  if (string) format(out, format) else out
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

  meta <- sprintf("v%s [%s]", archs4_version(x), archs4_creation_date(x))
  out <- paste(
    "=======================================================================\n",
    sprintf("%s (%s)\n", name(x), meta),
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
