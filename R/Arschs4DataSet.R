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

  out$nstudies <- length(unique(out$samples$dataset))
  out$sample_source <- attr(out$samples, "source")

  class(out) <- c("FacileArchs4DataSet", "FacileDataStore", class(out))
  out$samples <- FacileData::set_fds(out$samples, out)
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
  out <- paste(
    "=======================================================================\n",
    name(x), "\n",
    "-----------------------------------------------------------------------\n",
    sprintf("  path: %s\n", x$h5),
    sprintf("  studies: %s\n", prettyNum(x$nstudies, big.mark = ",")),
    sprintf("  samples: %s\n", prettyNum(nrow(x$samples), big.mark = ",")),
    sprintf("  genes: %s\n", prettyNum(nrow(x$features), big.mark = ",")),
    "  .....................................................................\n",
    sprintf("  sample source: %s\n", x$sample_source),
    "=======================================================================\n",
    sep = ""
  )
}

sample_search <- function(
  x,
  regex,
  meta_fields = NULL,
  scrub_search = FALSE,
  ignore.case = TRUE,
  remove_sc = x$remove_sc,
  threshold_sc = x$threshold_sc
) {
  assert_class(x, "FacileArchs4DataSet")
  out <- fuzzy_search(
    x$samples,
    regex = regex,
    meta_fields = meta_fields,
    scrub_search = scrub_search,
    ignore.case = ignore.case,
    remove_sc = remove_sc,
    threshold_sc = threshold_sc
  )
  out
}
