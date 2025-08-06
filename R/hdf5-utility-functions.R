# Utility functions to read archs4 specific files ------------------------------
.default_meta_fields <- function() {
  c(
    # "geo_accession", "series_id",
    "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")
}

.load_archs4_features <- function(path, species, verbose = FALSE, ...) {
  tictoc::tic("load features")
  out <- .hdf5_group_load_table(path, "meta/genes")
  tictoc::toc(quiet = !verbose)
  species <- match.arg(species, c("mouse", "human"))
  # in v2.2, the feature table has a few different columns, but there is a
  # subset of columns that are consistent
  dplyr::transmute(
    out,
    h5idx,
    # feature_id = ensembl_gene_id, # v2.2
    feature_id = ensembl_gene, # v2.2
    name = symbol,
    meta = biotype)
}

.load_archs4_samples <- function(
    path,
    columns = NULL,
    use_cache = TRUE,
    verbose = FALSE,
    ...
) {
  if (isTRUE(columns == "all")) {
    columns <- NULL
  } else {
    if (is.null(columns) || length(columns) == 0L) {
      columns <- c(
        "series_id",     # GSE160572
        "geo_accession", # GSM5899130
        "sample",        # GSM5899130
        "source_name_ch1",
        "title",
        "singlecellprobability",
        "alignedreads",
        "molecule_ch1", # cyotplasmic RNA, nuclear RNA, polyA RNA, etc.
        "library_source",    # transcriptomic, transcriptomic singlecell (not reliable)
        "library_selection", # CAGE, cDNA (vast majority), RACE, other
        "platform_id", # GPL24676, GPL20301, ...
        "instrument_model", # Illumina HiSeq 4000, etc.
        "data_processing",
        "extract_protocol_ch1",
        "characteristics_ch1"
      )
    }
  }
  tictoc::tic("sample load")
  out <- NULL
  source <- "h5"
  if (use_cache) {
    cfn <- archs4_cache_fn(path)$samples
    if (!file.exists(cfn)) {
      warning(
        "Samples cache not generated, consider running ",
        "`archs4_create_cache_files(", path, ")`"
      )
    } else {
      out <- arrow::read_parquet(cfn)
      bad.cols <- setdiff(columns, colnames(out))
      if (length(bad.cols) > 0) {
        warning(
          "This columns not present in samples [ignoring]: ",
          paste(bad.cols, collapse = ",")
        )
      }
      columns <- unique(c("h5idx", columns))
      out <- dplyr::select(out, dplyr::any_of(columns))
      source <- cfn
    }
  }

  if (is.null(out)) {
    out <- .hdf5_group_load_table(path, "meta/samples", columns = columns)
  }
  tictoc::toc(quiet = !verbose)
  attr(out, "source") <- source
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

# Even Lower Level Utility Functions -------------------------------------------
#' Read data from an HDF5 file "with caution"
#'
#' @description
#' This is a simple wrapper to [rhdf5::h5read()] which returns a default value
#' if the "data element" (specified by `name`) does not exists within `file`.
#' We use this file to read from ARCHS4 hdf5 files when we want to provide a
#' little insurance to the evolving nature of their data formats.
#'
#' For instance this function is used when we try to read something like
#' `"meta/reads_aligned"` because this information was not provided in earlier
#' versions of these datasets, however `"meta/genes"` may use [rhdf5::h5read()]
#' directly because this has been around since "the beginning"
#'
#' @importFrom rhdf5 h5read
#' @param default_value the default value to return if the entry isn't found,
#'   (Default: `NA`).
#' @param default_dim number of dimensions the result is expected to have,
#'   (Default: `1`).
.h5read <- function(file, name, index = NULL, start = NULL, stride = NULL,
                    block = NULL, count = NULL, compoundAsDataFrame = TRUE,
                    callGeneric = TRUE, read.attributes = FALSE, drop = FALSE,
                    ..., native = FALSE, s3 = FALSE, s3credentials = NULL,
                    default_value = NA, default_dim = 1L) {
  out <- try({
    rhdf5::h5read(file = file, name = name, index = index, start = start,
                  stride = stride, block = block, count = count,
                  compoundAsDataFrame = compoundAsDataFrame,
                  callGeneric = callGeneric, read.attributes = read.attributes,
                  drop = drop, ..., native = native, s3 = s3,
                  s3credentials = s3credentials)
  }, silent = TRUE)
  if (is(out, "try-error")) {
    ndim <- length(default_dim)
    if (ndim == 1L) {
      out <- rep(default_value, default_dim)
    } else {
      out <- array(default_value, default_dim)
    }
  }

  if (is.null(dim(out)) || length(dim(out)) == 1L) {
    na.it <- out %in% c("na", "NA", "null", "NULL")
    out[na.it] <- NA
  }

  out
}
