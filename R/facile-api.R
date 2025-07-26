#' @noRd
#' @export
default_assay.Archs4Client <- function(x, ...) {
  "counts"
}

#' Retrieve assay data
#' @export
#'
fetch_assay_data.archs4_client_facile_frame <- function(
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
