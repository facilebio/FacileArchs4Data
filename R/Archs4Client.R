#' Thin wrapper to the hdf5 files of an Archs4 dataset.
#'
#' Create a wrapper to the hdf5 files that are directly downloaded from the
#' archs4 website.
#'
#' To make your life easy, store the "equally versioned" mouse and human
#' datasets in a single directory with the value of `species` in the filename.
#' You can put new/old h5 files in different directories, then load the dataset
#' by specifying simply the species. See **Examples** for examples.
#'
#' @export
#' @param species `"human"` or `"mouse"`
Archs4Client <- R6::R6Class(
  "Archs4Client",
  public = list(
    #' @field path the filepath to the species-specific hdf5 data object
    path = NULL,

    #' @field samples the facile_frame of the samples in this object
    samples = NULL,

    #' @field features the tibble of features (genes) in this dataset.
    features = NULL,

    #' @field species is this thing mouse or human?
    species = NULL,

    #' @field remove_sc default value for the `remove_sc` parameters
    remove_sc = NULL,

    #' @description
    #' Create a new Archs4Client object
    #' @param species `"mouse"` or `"human"`
    #' @param ... stuff
    #' @param path the full path to the h5 file to load
    #' @param directory the parent directory that stores the species-specific
    #'   expression datasets. Defaults to `getOption("archs4.data_dir")`. To make
    #'   your life easier, you can set this option in your `~/.Rprofile`
    #' @param sample_columns wut?
    #' @param remove_sc default value for the `remove_sc` parameters
    #' @return a new Archs4Client
    #' @examples
    #' \dontrun{
    #' options(archs4.data_dir = "/path/to/directory/with/hdf5-count-files")
    #' a4m <- Archs4Client$new("mouse")
    #' a4h <- Archs4Client$new("/path/to/human_gene_v2.2.h5")
    #' }
    initialize = function(species, ..., path = NULL,
                          directory = getOption("archs4.data_dir", NULL),
                          remove_sc = TRUE,
                          sample_columns = NULL) {
      self$remove_sc <- assert_flag(remove_sc)
      if (is.null(path)) {
        assert_directory_exists(directory, "r")
        assert_string(species)
        path <- dir(directory, species, full.names = TRUE)
      }
      assert_string(path)
      path <- normalizePath(assert_file_exists(path, "r", extension = "h5"))
      if (grepl("mouse", basename(path), ignore.case = TRUE)) {
        self$species <- "mouse"
      } else {
        self$species <- "human"
      }
      self$path <- path

      use_cache <- getOption("archs4.use_cache", TRUE)
      cached <- .pkgcache[[path]]
      in_cache <-
        is.list(cached) &&
        all(c("samples", "features") %in% names(cached))

      if (!in_cache) {
        cached <- list(
          samples = .load_archs4_samples(path, columns = sample_columns),
          features = .load_archs4_features(path, self$species))
        if (use_cache) {
          .pkgcache[[path]] <- cached
        }
      }

      self$samples <- set_fds(cached$samples, self)
      self$features <- cached$features
      self
    },

    #' @description
    #' Find samples that match a regular expression in any of the descriptive
    #' columns.
    #'
    #' This function is fashioned after the archs4py.meta() function.
    #'
    #' @seealso https://github.com/MaayanLab/archs4py/blob/main/archs4py/data.py#L27
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
    search = function(regex, meta_fields = NULL,
                      scrub_search = FALSE, ignore.case = TRUE,
                      remove_sc = self$remove_sc,
                      threshold_sc = 0.5) {
      out <- fuzzy_search(
        self$samples,
        regex,
        meta_fields = meta_fields,
        scrub_search = scrub_search,
        ignore.case = ignore.case,
        remove_sc = self$remove_sc,
        threshold_sc = 0.5
      )
      out
    },

    #' @description
    #' Return a table of studies (series_id) in here, optionally dropping
    #' those that look singlecell-like
    #'
    #' @param threshold_sc the probability threshold to indicate wheter sample
    #'   is from a singlecell dataset, default: 0.5
    #' @param threshold_sc_likely the proportion of samples in a study that
    #'   must look like singlecell data in order to classify the study as
    #'   singlecell. default: 0.2
    #' @param sample_search If a search string is provided, this will trigger
    #'   a sample search via the `$search()` method. The study table returned
    #'   will be augmented with a `search_hits` column that tells you how many
    #'   samples within the study were a hit for your search.
    #' @param scrub_search,meta_fields,ignore.case parameters passed to `search` function.
    #' @return a tibble with 1 row per geo series (GSE) / study with the
    #'   following columns:
    #'   * `nsamples`: number of samples in the study
    #'   * `nsc`: number of samples that look like single cell, given the
    #'     `threshold_sc` parameter
    #'   * `sc_likely`: if the number of samples in the study that look like
    #'     they are singlecell exceeds `threshold_sc_likely`, this is set to
    #'     `TRUE`
    #'   * `search_hits`: The number of samples in the series that were a hit
    #'     for the `sample_search` query. If no `sample_search` was provided,
    #'     this is just `0`.
    studies = function(
        sample_search = NULL,
        meta_fields = .default_meta_fields(),
        scrub_search = FALSE,
        ignore.case = TRUE,
        threshold_sc = 0.5,
        threshold_sc_likely = 0.2
      ) {
      assert_number(threshold_sc, lower = 0.01, upper = 1)

      hits <- NULL
      if (checkmate::test_string(sample_search)) {
        hits <- self$search(sample_search, meta_fields, ignore.case) |>
          dplyr::summarize(search_hits = dplyr::n(), .by = "series_id")
      }
      smry <- self$samples |>
        dplyr::summarize(
          nsamples = dplyr::n(),
          nsc = sum(singlecellprobability >= threshold_sc),
          likely_sc = nsc / nsamples >= 0.2,
          search_hits = 0L,
          .by = "series_id"
        )
      if (!is.null(hits) && nrow(hits) > 0L) {
        smry <- smry |>
          dplyr::select(-search_hits) |>
          dplyr::left_join(hits, by = "series_id") |>
          dplyr::mutate(search_hits = ifelse(is.na(search_hits), 0L, search_hits))
      }
      smry
    }
  ),
  private = list(
    load_samples = function(path) {}
  )
)

