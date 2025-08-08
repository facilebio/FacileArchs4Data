.infer_species_from_filename <- function(x) {
  assert_character(x)
  is_mouse <- grepl("mouse", basename(x))
  is_human <- grepl("human", basename(x))
  ok <- xor(is_mouse, is_human)
  if (!all(ok)) {
    stop("Couldn't infer if file is human or mouse: ", basename(x)[!ok])
  }
  ifelse(is_mouse, "mouse", "human")
}

#' List the archs4 files in a directory
#'
#' archs4-related files need to have `"moue"` or `"human"` in their
#' filename.
#'
#' @export
#' @param x the directory to look inside
archs4_dir_info <- function(x = getOption("archs4.data_dir"), ...) {
  assert_directory_exists(x, "r")
  regex.h5 <- "^\\d{4}-\\d{2}-\\d{2}.*\\.(h5|hdf5)$"

  # This pipeline doesn't break if no files are found
  dinfo.h5 <- fs::dir_info(x) |>
    dplyr::mutate(fn = basename(path)) |>
    dplyr::filter(grepl(regex.h5, fn)) |>
    dplyr::filter(type == "file") |>
    dplyr::select(fn, path, birth_time) |>
    dplyr::mutate(
      species = .infer_species_from_filename(fn),
      arcsh4_type = "release",
      .after = "fn"
    )
  if (nrow(dinfo.h5) == 0L) return(dinfo.h5)

  meta.h5 <- archs4_meta(dinfo.h5$path, c("version", "creation_date"))

  h5.info <- dinfo.h5 |>
    dplyr::bind_cols(meta.h5) |>
    dplyr::mutate(
      sample_cache_path = sapply(path, function(y) {
        archs4_cache_fn(y, outdir = x)$sample
      }),
      sample_cache_exists = file.exists(sample_cache_path)
    ) |>
    dplyr::relocate(dplyr::ends_with("path"), .after = dplyr::last_col()) |>
    dplyr::arrange(dplyr::desc(birth_time))
  h5.info
}
