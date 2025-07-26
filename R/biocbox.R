#' @noRd
#' @export
biocbox.archs4_client_facile_frame <- function(
    x, class = "list", assay_name = "counts", features = NULL,
    sample_covariates = NULL, with_h5idx = FALSE, ...) {
  a4c <- assert_class(fds(x), "Archs4Client")
  class <- assert_choice(class, c("list", "DGEList", "SummarizedExperiment"))
  assert_integerish(x$h5idx)
  if (is.null(features)) {
    features <- a4c$features
  }
  assert_tibble(features)
  assert_subset(c("h5idx", "feature_id", "meta", "name"),  names(features))

  counts <- rhdf5::h5read(
    a4c$path,
    "data/expression",
    index = list(x$h5idx, features$h5idx)) |>
    t()
  rownames(counts) <- features$feature_id

  if (!any(c("sample_id", "sample") %in% colnames(x))) {
    warning(
      "neither `sample_id` or `sample` column found in sample frame, ",
      "creating pseudo id's")
    x$sample_id <- sprintf("sample_%d", seq_len(ncol(counts)))
  }
  if ("sample_id" %in% colnames(x)) {
    colnames(counts) <- x$sample_id
  } else {
    colnames(counts) <- x$sample
  }

  # Let's work with the sample frame we have, not the one we think we have, ie.
  # we may have manipulated / cleaned this up before passing it in here.
  if (!"dataset" %in% names(x) && "series_id" %in% names(x)) {
    x <- dplyr::rename(x, dataset = series_id)
  }
  if (!with_h5idx) {
    x$h5idx <- NULL
  }

  out <- list(
    counts = counts,
    features = features |>
      as.data.frame() |>
      dplyr::mutate(
        feature_type = "ensgid",
        source = "ensembl_v107"),
    samples = as.data.frame(x))
  rownames(out$features) <- rownames(counts)
  rownames(out$samples) <- colnames(counts)

  if (class == "DGEList") {
    reqpkg("edgeR")
    out <- edgeR::DGEList(
      counts = out$counts,
      genes = out$features,
      samples = out$samples) |>
      edgeR::calcNormFactors()
  } else if (class == "SummarizedExperiment") {
    reqpkg("SummarizedExperiment")
    out <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = out$counts),
      rowData = out$features,
      colData = out$samples)
  }

  out
}
