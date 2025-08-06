# Implementing same functions required in FacileBiocData to make some magic


#' Shows tibble of assay information
#' @noRd
#' @export
assay_info.Archs4Client <- function(x, assay_name = NULL, ...) {
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


