# Implementing same functions required in FacileBiocData to make some magic


#' Shows tibble of assay information
#' @noRd
#' @export
assay_info.FacileArchs4DataSet <- function(x, assay_name = NULL, ...) {
  info <- dplyr::tribble(
    # do we want to add TPM here? -- we can technically calculate on fly
    ~assay,   ~assay_type, ~feature_type, ~description,         ~nfeatures,       ~storatge_mode,
    "counts", "rnaseq",    "ensgid",      "counts from archs4", nrow(x$features), "integer"
  )
  info
}


#' @noRd
#' @export
default_assay.Archs4Client <- function(x, ...) {
  "counts"
}

#' @noRd
#' @export
samples.FacileArchs4DataSet <- function(x, ..., drop_covariates = TRUE) {
  assert_class(x, "FacileArchs4DataSet")
  out <- x$samples
  if (drop_covariates) {
    out <- dplyr::select(out, dataset, sample_id)
  }
  out
}

#' Retrieve assay data
#' @export
#'
fetch_assay_data.archs4_facile_frame <- function(
  x,
  features = NULL,
  samples = NULL,
  assay_name = NULL,
  normalized = FALSE,
  batch = NULL,
  main = NULL,
  as.matrix = FALSE,
  drop_samples = TRUE,
  ...,
  prior.count = 3,
  subset.threshold = 700,
  aggregate = FALSE,
  aggregate.by = "ewm",
  verbose = FALSE
) {
}

#' Retrieve (wide) assay data
#' @export
with_assay_data.archs4_client_facile_frame <- function(
  x,
  features,
  assay_name = NULL,
  normalized = TRUE,
  aggregate = FALSE,
  aggregate.by = "ewm",
  spread = TRUE,
  with_assay_name = FALSE,
  ...,
  verbose = FALSE,
  .fds = fds(x)
) {
}
