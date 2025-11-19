# Utility functions to read archs4 specific files ------------------------------

#' @noRd
.default_meta_fields <- function() {
  c(
    # "geo_accession", "series_id",
    "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")
}

#' @noRd
.load_archs4_metadata <- function(x, ..., as_list = TRUE) {
  meta <- .hdf5_group_load_table(x, "meta/info")
  meta[["h5idx"]] <- NULL
  names(meta) <- janitor::make_clean_names(names(meta))
  meta$creation_date <- as.Date(meta$creation_date, "%m/%d/%Y")
  if (as_list) as.list(meta) else meta
}

#' Guess the speices from the hdf5 file ensembl identifiers used
#' @noRd
.hdf5_guess_species <- function(path, ...) {
  checkmate::assert_file_exists(path)
  sapply(path, function(p) {
    xx <- .hdf5_group_load_table(p, "meta/genes", rows = 1)
    out <- NULL
    if (grepl("ENSG\\d+$", xx$ensembl_gene)) {
      out <- "human"
    } else if (grepl("ENSMUSG\\d+$", xx$ensembl_gene)) {
      out <- "mouse"
    } else {
      stop("Can't figure out species")
    }
    out
  })
}

#' @noRd
.load_archs4_features <- function(path, species = NULL, verbose = FALSE, ...) {
  if (is.null(species)) {
    species <- .hdf5_guess_species(path)
  }
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

#' @noRd
.load_archs4_samples <- function(
    path,
    columns = NULL,
    rewrite_cache = FALSE,
    sample_cache_path = NULL,
    verbose = FALSE,
    cache_dir = NULL,
    ...,
    .duckdb = FALSE
) {
  assert_file_exists(path, "r")
  assert_logical(rewrite_cache)
  assert_logical(.duckdb)

  req.columns <- c(
    "h5idx",
    "dataset",
    "sample_id",
    "alignedreads",
    "singlecellprobability",
    "characteristics_ch1",
    "source_name_ch1",
    "title"
  )

  if (isTRUE(columns == "all")) {
    columns <- NULL
  } else {
    if (is.null(columns) || length(columns) == 0L) {
      columns <- c(
        req.columns,
        "molecule_ch1", # cyotplasmic RNA, nuclear RNA, polyA RNA, etc.
        "extract_protocol_ch1",
        "library_source",    # transcriptomic, transcriptomic singlecell (not reliable)
        "library_selection", # CAGE, cDNA (vast majority), RACE, other
        "platform_id", # GPL24676, GPL20301, ...
        "instrument_model", # Illumina HiSeq 4000, etc.
        "data_processing"
      )
    }
  }
  columns <- unique(c(req.columns, columns))

  tictoc::tic("sample load")
  out <- NULL

  source <- "h5"

  if (is.null(sample_cache_path)) {
    sample_cache_path <- archs4_cache_fn(path, cache_dir = cache_dir)$samples
  }

  cdir <- dirname(sample_cache_path)
  if (!dir.exists(cdir)) {
    if (!dir.create(cdir)) {
      stop("Could not make cache directory: ", cdir)
    }
  }

  if (!file.exists(sample_cache_path) || rewrite_cache) {
    out <- .hdf5_group_load_table(path, "meta/samples", columns = NULL) |>
      dplyr::rename(dataset = series_id, sample_id = geo_accession) |>
      dplyr::relocate(out, dplyr::all_of(req.columns), .before = 1L)
    arrow::write_parquet(out, sample_cache_path)
  } else {
    source <- sample_cache_path
    if (.duckdb) {
      out <- duckdb::duckdb() |>
        DBI::dbConnect() |>
        dplyr::tbl(sprintf("read_parquet('%s')", sample_cache_path))
    } else {
      out <- arrow::read_parquet(sample_cache_path)
    }
  }

  if (!inherits(out, "tbl_dbi") && !is.null(columns)) {
    bad.cols <- setdiff(columns, colnames(out))
    if (length(bad.cols) > 0) {
      warning(
        "This columns not present in samples [ignoring]: ",
        paste(bad.cols, collapse = ",")
      )
    }
    out <- dplyr::select(out, dplyr::any_of(columns))
  }
  tictoc::toc(quiet = !verbose)

  attr(out, "source") <- source
  out
}

#' @noRd
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

#' @noRd
.hdf5_group_load_table <- function(
  x,
  group = "meta/genes",
  columns = NULL,
  rows = NULL
) {
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

  index <- if (checkmate::test_integerish(rows)) as.integer(rows) else NULL

  info <- lapply(keep, function(vname) {
    info <- dplyr::filter(group_info, .data$name == .env$vname)
    vals <- rhdf5::h5read(
      x,
      paste0(group, "/", info$name),
      index = index,
      bit64conversion = "double"
    )
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
