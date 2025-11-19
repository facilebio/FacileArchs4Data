#' List the archs4 files in a directory
#'
#' archs4-related files need to have `"moue"` or `"human"` in their
#' filename.
#'
#' Note that the `creation_date` column is extracted from the hdf5 file, so
#' this is the source of truth for the version of this data file. It is
#' not the same as the date listed on the website, for instance on 2025-11-14
#' the release of the gene-level datafile is 2025-11-11, but the hdf5 file
#' has a creation date of 2025-11-06
#'
#' @export
#' @param x the directory to look inside
archs4_dir_info <- function(
  x = getOption("archs4.data_dir"),
  cache_dir = file.path(x, "cache"),
  hdf5_suffix = "h5",
  ...
) {
  assert_directory_exists(x, "r")
  # regex.h5 <- "^\\d{4}-\\d{2}-\\d{2}.*\\.(h5|hdf5)$"
  regex.h5 <- sprintf(".*\\.(%s)", paste(hdf5_suffix, collapse = "|"))

  # This pipeline doesn't break if no files are found
  dinfo.h5 <- fs::dir_info(x) |>
    dplyr::mutate(fn = basename(path)) |>
    dplyr::filter(grepl(regex.h5, fn)) |>
    dplyr::filter(type == "file") |>
    dplyr::select(fn, path) |>
    dplyr::mutate(
      # species = .infer_species_from_filename(fn),
      species = .hdf5_guess_species(path),
      archs4_type = "release",
      .after = "fn"
    )

  if (nrow(dinfo.h5) == 0L) {
    return(dinfo.h5)
  }

  meta.h5 <- archs4_meta(dinfo.h5$path, c("version", "creation_date"))

  h5.info <- dinfo.h5 |>
    dplyr::bind_cols(meta.h5) |>
    dplyr::relocate(
      dplyr::all_of(colnames(meta.h5)),
      .before = "archs4_type"
    ) |>
    dplyr::mutate(
      sample_cache_path = sapply(path, function(y) {
        archs4_cache_fn(y, cache_dir = cache_dir)$sample
      }),
      sample_cache_exists = file.exists(sample_cache_path)
    ) |>
    dplyr::relocate(dplyr::ends_with("path"), .after = dplyr::last_col()) |>
    dplyr::arrange(dplyr::desc(creation_date))

  h5.info
}

#' @noRd
.infer_species_from_filename <- function(x) {
  .Deprecated(".hdf5_guess_species")
  assert_character(x)
  is_mouse <- grepl("mouse", basename(x))
  is_human <- grepl("human", basename(x))
  ok <- xor(is_mouse, is_human)
  if (!all(ok)) {
    stop("Couldn't infer if file is human or mouse: ", basename(x)[!ok])
  }
  ifelse(is_mouse, "mouse", "human")
}
