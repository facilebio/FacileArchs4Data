#' Thin wrapper to the hdf4 files of an Archs4 dataset.
#'
#' To make your life easy, store the "equally versioned" mouse and human
#' datasets in a single directory with the value of `species` in the filename.
#' You can put new/old h5 files in different directories, then load the dataset
#' by specifying simply the species. See **Examples** for examples.
#'
#' @export
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

    #' @description
    #' Create a new Archs4Client object
    #' @param species `"mouse"` or `"human"`
    #' @param ... stuff
    #' @param path the full path to the h5 file to load
    #' @param directory the parent directory that stores the species-specific
    #'   expression datasets. Defaults to `getOption("archs4.data_dir")`. To make
    #'   your life easier, you can set this option in your `~/.Rprofile`
    #' @param sample_columns wut?
    #' @return a new Archs4Client
    #' @examples
    #' \dontrun{
    #' options(archs4.data_dir = "/path/to/directory/with/hdf5-count-files")
    #' a4m <- Archs4Client$new("mouse")
    #' a4h <- Archs4Client$new("/path/to/human_gene_v2.2.h5")
    #' }
    initialize = function(species = c("human", "mouse"), ..., path = NULL,
                          directory = getOption("archs4.data_dir", NULL),
                          sample_columns = NULL) {
      if (is.null(path)) {
        assert_directory_exists(directory, "r")
        species <- match.arg(species)
        path <- dir(directory, species, full.names = TRUE)
      }
      checkmate::assert_string(path)
      checkmate::assert_file_exists(path, "r", extension = "h5")
      self$species <- if (grepl("mouse", basename(path), ignore.case = TRUE)) {
        "mouse"
      } else {
        "human"
      }
      self$path <- path
      self$samples <- set_fds(
        .load_archs4_samples(path, columns = sample_columns),
        self)
      self$features <- .load_archs4_features(path, self$species)

    }
  ),
  private = list(
    load_samples = function(path) {}
  )
)

.load_archs4_features <- function(path, species) {
  out <- .hdf5_group_load_table(path, "meta/genes")
  species <- match.arg(species, c("mouse", "human"))
  # in v2.2, the feature table has a few different columns, but there is a
  # subset of columns that are consistent
  dplyr::transmute(
    out,
    h5idx,
    feature_id = ensembl_gene_id,
    name = symbol,
    meta = biotype)
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
  out <- .hdf5_group_load_table(path, "meta/samples", columns = columns)
  class(out) <- c("archs4_client_facile_frame", class(out))
  out
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
