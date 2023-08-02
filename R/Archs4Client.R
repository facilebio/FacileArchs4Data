#' Thin wrapper to Archs4 dataset
#'
#' @export
Archs4Client <- R6::R6Class(
  "Archs4Client",
  public = list(
    path = NULL,
    samples = NULL,
    features = NULL,

    initialize = function(path, ..., sample_columns = NULL) {
      checkmate::assert_file_exists(path, "r", extension = "h5")
      self$path <- path
      self$samples <- .load_archs4_samples(path, columns = sample_columns)
      self$features <- .load_archs4_features(path)
    }
  ),
  private = list(
    load_samples = function(path) {}
  )
)

.load_archs4_features <- function(path) {
  .hdf5_group_load_table(path, "meta/genes")
}

.load_archs4_samples <- function(path, columns = NULL) {
  if (is.null(columns) || length(columns) == 0L) {
    columns = c(
      "series_id",     # GSE160572
      "geo_accession", # GSM5899130
      "sample",        # GSM5899130
      "source_name_ch1",
      "title",
      "singlecellprobability",
      "readsaligned",
      "molecule_ch1", # cyotplasmic RNA, nuclear RNA, polyA RNA, etc.
      "library_source",    # transcriptomic, transcriptomic singlecell (not reliable)
      "library_selection", # CAGE, cDNA (vast majority), RACE, other
      "platform_id", # GPL24676, GPL20301, ...
      "instrument_model", # Illumina HiSeq 4000, etc.
      "data_processing",
      "extract_protocol_ch1",
      "characteristics_ch1")
  }
  .hdf5_group_load_table(path, "meta/samples", columns = columns)
}

.hdf5_group_info <- function(x, group = "meta/genes") {
  if (checkmate::test_string(x)) {
    checkmate::assert_file_exists(x, "r", extension = "h5")
    x <- rhdf5::H5Fopen(x)
    on.exit(rhdf5::H5Fclose(x), add = TRUE)
  }
  checkmate::assert_class(x, "H5IdComponent")
  hdf5_group <- rhdf5::H5Gopen(x, group)
  on.exit(rhdf5::H5Gclose(hdf5_group), add = TRUE)
  group_info <- rhdf5::h5ls(hdf5_group)
  group_info
}

.hdf5_group_load_table <- function(x, group = "meta/genes", columns = NULL) {
  group_info <- .hdf5_group_info(x, group)
  if (is.character(columns)) {
    invalid <- setdiff(columns, group_info$name)
    if (length(invalid)) {
      warning("Unknown columns: ", paste(invalid, collapse = ","))
    }
    keep <- intersect(columns, group_info$name)
    if (length(keep) == 0L) stop("No columns will be returned")
    group_info <- dplyr::filter(group_info, .data$name %in% .env$keep)
  } else {
    keep <- group_info$name
  }

  info <- lapply(keep, function(vname) {
    info <- dplyr::filter(group_info, .data$name == .env$vname)
    vals <- rhdf5::h5read(x, paste0(group, "/", info$name),
                          bit64conversion = "double")
    if (info$dclass == "INTEGER") {
      vals <- as.integer(vals)
    } else if (info$dclass == "FLOAT") {
      vals <- as.numeric(vals)
    } else {
      vals <- as.character(vals)
    }
    vals
  })
  names(info) <- keep
  out <- dplyr::as_tibble(info)
  dplyr::mutate(out, h5idx = seq(nrow(out)), .before = 1L)
}
