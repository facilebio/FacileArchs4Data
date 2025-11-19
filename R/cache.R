#' Create sample cache for hdf5 files
#'
#' You can point to a specific file, or the arcsh4 directory. In the latter
#' case a cache file will be generated for all hdf5 files that do not have one.
#'
#' @export
#' @param hdf5 explicit path to the hdf5 file to create cache for
#' @param datadir
archs4_create_sample_cache <- function(
  hdf5 = NULL,
  datadir = getOption("archs4.data_dir"),
  cache_dir = NULL,
  ...,
  create_cache_dir = TRUE
) {
  if (is.null(hdf5)) {
    checkmate::assert_directory_exists(datadir)
    ainfo <- archs4_dir_info()
  }
}


# Original Cache Realted Functions ---------------------------------------------

#' File paths for cached objects created from an archs4 hdf5 file.
#'
#' The name of the cache file that is generated does not relate to the file
#' name of the hdf5 file, but rather the creation date and species information
#' that is stored from the metadata within the hdf5 file.
#'
#' @export
#' @param x the path to an archs4 data file
#' @param outdir the directory where the cache files are expected to be.
#'   Defaults to `dirname(x)`
#' @examples
#' archs4_cache_fn("~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5")
archs4_cache_fn <- function(x, cache_dir = NULL, ...) {
  checkmate::assert_string(x) # one hdf5 file at a time
  checkmate::assert_file_exists(x, extension = "h5")
  if (is.null(cache_dir)) {
    cache_dir <- file.path(dirname(x), "cache")
  }

  species <- .hdf5_guess_species(x)
  creation_date <- archs4_creation_date(x)

  bname <- sub("\\.(h5|hdf5)$", "", basename(x))
  outfn <- list(
    # samples = sprintf("%s__%s-samples.parquet", creation_date, species)
    samples = paste0(bname, ".sample_cache.parquet")
  )

  lapply(outfn, function(fn) file.path(cache_dir, fn))
}

#' Save the samples and the gene feature info for an archs4 dataset as paruqet
#' files to make loading quicker.
#'
#' Loading the samples table from the hdf5 takes ~25s vs 1s from parquet
#'
#' @param x the path to the hdf5 file
#'
#' @examples
#' h5fn <- "~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5"
#' archs4_create_cache_files(h5fn)
#' fns <- archs4_cache_fn(h5fn)
#' sapply(fns, file.exists)
archs4_create_cache_files <- function(
  x = getOption("archs4.data_dir"),
  cache_dir = dirname(x),
  ...
) {
  if (FALSE) {
    x <- getOption("archs4.data_dir")
    cache_dir <- dirname(x)
  }
  checkmate::assert_string(x)
  stopifnot(
    "input is file or directory" = file.exists(x)
  )
  if (checkmate::test_directory(x)) {
    inputs <- archs4_dir_info(x)
  } else {
    checkmate::assert_file_exists(x, extension = c("h5", "hdf5"))
    inputs <- archs4_dir_info(dirname(x), hdf5_suffix = fs::path_ext(x)) |>
      dplyr::filter(fn %in% basename(x))
  }
  checkmate::assert_file_exists(x, extension = "h5")
  checkmate::assert_directory_exists(outdir, "w")
  if (!exists(cache_dir)) {
    dir.create(cache_dir)
  }
  outfn <- archs4_cache_fn(x, cache_dir = cache_dir)

  # samples
  # xsamples <- .hdf5_group_load_table(x, "meta/samples", columns = NULL)
  xsamples <- .load_archs4_samples(x, cache_dir = cache_dir)
  arrow::write_parquet(xsamples, outfn$samples)
  invisible(outfn)
}

#' Pull down gse experiment information
#' @export
archs4_cache_experiment_details <- function(
  x,
  ncores = parallel::detectCores() - 1,
  outdir = NULL,
  ...
) {
  checkmate::assert_class(x, "FacileArchs4DataSet")
  if (is.null(outdir)) {
    outdir <- sub("h5", "cache/series_info", x$h5)
  }
  if (!dir.exists(outdir)) {
    stopifnot(dir.create(outdir, recursive = TRUE))
  }
  checkmate::assert_directory(outdir, "w")
  gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

  geoids <- samples(x) |>
    dplyr::count(dataset) |>
    dplyr::collect() |>
    tidyr::separate_longer_delim(dataset, ",")

  # geoids <- head(geoids, 20)
  # geoids$dl <- parallel::mclapply(head(geoids$dataset, 20), function(gid) {
  geoids$dl <- lapply(geoids$dataset, function(gid) {
    tryCatch(
      {
        GEOquery::getGEOfile(
          gid,
          destdir = outdir,
          AnnotGPL = FALSE,
          amount = "quick"
        )
      },
      error = function(e) NA_character_
    )
  })
  # }, mc.preschedule = TRUE, mc.cores = 1)
  geoids
}
