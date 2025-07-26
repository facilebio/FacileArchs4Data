# Some cache-realted functions.


#' File paths for cached objects created from an archs4 hdf5 file
#'
#' @export
#' @param x the path to an archs4 data file
#' @param outdir the directory where the cache files are expected to be.
#'   Defaults to `dirname(x)`
#' @examples
#' archs4_cache_fn("~/workspace/data/archs4/v2.6/mouse_gene_v2.6.h5")
archs4_cache_fn <- function(x, outdir = dirname(x), ...) {
  checkmate::assert_file_exists(x, extension = "h5")
  checkmate::assert_directory_exists(outdir, "w")
  outfn <- list(
    samples = sub("h5$", "samples.parquet", basename(x))
  )
  lapply(outfn, function(fn) file.path(outdir, fn))
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
