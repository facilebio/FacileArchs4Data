#' The `characteristics_ch1` column often has key:value attributes worth parsing.
#'
#' Often times the archs4 `characteristics_ch1` sample column looks like
#' `"cell line: K562,cell type: Immortalized epithelial cells,genotype: WT, ...`
#' Let's bust this up into usable data.frame columns.
#'
#' TODO: Fix this parsing to accomodate for multiple comma separated values in
#' on keyvalue par. For instance, the mouse GSE198071 experiment has
#' characteristics that looks like this:
#'
#' `"tissue: Tibialis anterior muscle,genotype: mfn1-/-, mfn2-/-,age: 6 weeks"`
#'
#' This gets into the appropriate columns (tissue, genotype, age), however the
#' genotype column only has mfn1-/- values. The mfn2-/- gets lost
#'
#' @export
#' @param x the samples archs4 data.frame
#' @param column the name of the column to key:value extract, default
#'   is `"characteristics_ch1"`
#' @param delim_var the delimeter used to split key:value pairs, default: `","`
#' @param delim_key the delimeter to split key/value pairs, default: `":"`
#' @param clean_names use [janitor::clean_names()] to clean up the key variable
#'   names expanded from `column`.
#' @param keep_column return the `x[[column]]` column in the output, defaults
#'   to `TRUE`
#' @param append if TRUE (default) returns `x` with the new columns, otherwise
#'   just returns the metadata
#' @param ... stuff and things
#' @return wide tibble of metadata variables from `x[[column]]`
#' @examples
#' a4h <- FacileArchs4Data::FacileArchs4DataSet("human")
#' info.raw <- samples(a4h) |> dplyr::filter(dataset == "GSE212337")
#' info <- archs4_split_metadata(info.raw)
#'
#' # This one has values w/ multiple ',' in genotype. I couldn't figure out
#' # how to encode the tidyr::split_* functions to accomodate for this, which
#' # is why the key<>value parsing is done inside loops
#' a4m <- FacileArchs4Data::FacileArchs4DataSet("mouse")
#' info.raw <- samples(a4m) |> dplyr::filter(dataset == "GSE198071")
#' info <- archs4_split_metadata(info.raw)
#' info |> select(h5idx:alignedreads)
#' table(info$genotype)
archs4_split_metadata <- function(
  x,
  column = "characteristics_ch1",
  delim_var = ",",
  delim_key = ":",
  append = TRUE,
  clean_names = TRUE,
  keep_column = !append,
  ...
) {
  checkmate::assert_character(x[["dataset"]])
  checkmate::assert_character(x[["sample_id"]])
  checkmate::assert_character(x[[column]])
  if (nrow(x) == 0L) return(x)

  if (FALSE) {
    info <- x |>
      dplyr::select(dataset, sample_id, dplyr::all_of(column)) |>
      tidyr::separate_longer_delim(column, delim = delim_var) |>
      tidyr::separate_wider_delim(
        column,
        delim = delim_key,
        names = c("key", "value"),
        too_many = "merge",
        too_few = "align_start"
      ) |>
      dplyr::mutate(
        key = janitor::make_clean_names(key, allow_dupes = TRUE),
        value = trimws(value)) |>
      tidyr::pivot_wider(names_from = "key", values_from = "value")
  }
  kv <- strsplit(x[[column]], delim_key)
  info <- lapply(seq_along(kv), function(i) {
    crap <- trimws(kv[[i]])
    if (length(crap) == 2L) {
      keys <- crap[1L]
      vals <- crap[2L]
    } else {
      n <- 2:(length(crap) - 1)
      keys <- unname(c(crap[1L], sub(".*,", "", crap[n])))
      vals <- sapply(strsplit(crap[-1L], delim_var), function(v) {
        if (length(v) == 1L) return(v)
        paste(trimws(head(v, -1L)), collapse = delim_var, sep = "")
      })
    }
    names(vals) <- janitor::make_clean_names(keys)
    dplyr::as_tibble(as.list(vals))
  }) |> dplyr::bind_rows()
  info <- dplyr::bind_cols(dplyr::select(x, dataset, sample_id), info)
  if (keep_column) {
    og <- dplyr::distinct(x, dataset, sample_id, .data[[column]])
    info <- dplyr::left_join(info, og, by = c("dataset", "sample_id"))
  }

  sample.count <- dplyr::count(info, dataset, sample_id)
  if (any(sample.count$n) > 1) {
    stop("split operation created duplicate rows for some samples")
  }

  if (append) {
    icols <- setdiff(colnames(info), c("dataset", "sample_id"))
    axe.cols <- intersect(icols, colnames(x))

    out <- x |>
      dplyr::select(-dplyr::all_of(axe.cols)) |>
      dplyr::left_join(info, by = c("dataset", "sample_id")) |>
      dplyr::relocate(dplyr::all_of(icols), .after = "sample_id")
    attr(out, "parsed") <- setdiff(icols, "characteristics_ch1")
  } else {
    out <- info
  }

  out
}
