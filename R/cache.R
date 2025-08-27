# Some cache-realted functions.

#' File paths for cached objects created from an archs4 hdf5 file
#'
#' @export
#' @param x the path to an archs4 data file
#' @param outdir the directory where the cache files are expected to be.
#'   Defaults to `dirname(x)`
#' @examples
#' archs4_cache_fn("~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5")
archs4_cache_fn <- function(x, outdir = NULL, ...) {
  checkmate::assert_string(x) # one hdf5 file at a time
  checkmate::assert_file_exists(x, extension = "h5")
  if (!is.null(outdir)) assert_directory_exists(outdir, "w")
  # outfn <- list(
  #   samples = sub("h5$", "samples.parquet", basename(x))
  # )
  species <- .infer_species_from_filename(x)
  basefn <- sub(".*?__", "", basename(x)) # strip off the data on the file
  outfn <- list(
    samples = sprintf(
      "%s__%s.samples.parquet",
      archs4_creation_date(x),
      sub("\\.(h5|hdf5)$", "", basefn)
    )
  )
  if (!is.null(outdir)) {
    outfn <- lapply(outfn, function(fn) file.path(outdir, fn))
  }
  outfn
}

#' Save the samples and the gene feature info for an archs4 dataset as paruqet
#' files to make load quicker.
#'
#' Loading the samples table from the hdf5 takes ~25s vs 1s from parquet
#' @param x the path to the hdf5 file
#' @examples
#' h5fn <- "~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5"
#' archs4_create_cache_files(h5fn)
#' fns <- archs4_cache_fn(h5fn)
#' sapply(fns, file.exists)
archs4_create_cache_files <- function(x, outdir = dirname(x), ...) {
  if (FALSE) {
    x <- "~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5"
    outdir <- dirname(x)
  }
  checkmate::assert_file_exists(x, extension = "h5")
  checkmate::assert_directory_exists(outdir, "w")
  outfn <- archs4_cache_fn(x, outdir = outdir)

  # samples
  xsamples <- .hdf5_group_load_table(x, "meta/samples", columns = NULL)
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
    tryCatch({
      GEOquery::getGEOfile(gid, destdir = outdir, AnnotGPL = FALSE, amount = "quick")
    }, error = function(e) NA_character_)
  })
  # }, mc.preschedule = TRUE, mc.cores = 1)
  geoids
}
