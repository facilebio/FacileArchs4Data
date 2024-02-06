#' @noRd
#' @export
biocbox.archs4_client_facile_frame <- function(
    x, class = "list", assay_name = "counts", features = NULL,
    sample_covariates = NULL, ...) {
  a4c <- assert_class(fds(x), "Archs4Client")
  class <- assert_choice(class, c("list", "DGEList", "SummarizedExperiment"))
  assert_integerish(x$h5idx)
  if (is.null(features)) {
    features <- a4c$features
  }
  assert_tibble(features)
  assert_subset(names(features), c("h5idx", "ensembl_gene_id", "biotype", "symbol"))

  counts <- rhdf5::h5read(
    a4c$path,
    "data/expression",
    index = list(x$h5idx, features$h5idx)) |>
    t()
  rownames(counts) <- features$ensembl_gene_id
  colnames(counts) <- x$sample

  # Let's work with the sample frame we have, not the one we think we have, ie.
  # we may have manipulated / cleaned this up before passing it in here.
  if (!"sample_id" %in% names(x) && "sample" %in% names(x)) {
    x <- dplyr::rename(x, sample_id = sample)
  }
  if (!"dataset" %in% names(x) && "series_id" %in% names(x)) {
    x <- dplyr::rename(x, dataset = series_id)
  }

  out <- list(
    counts = counts,
    features = features |>
      as.data.frame() |>
      dplyr::transmute(
        feature_id = ensembl_gene_id,
        feature_type = "ensgid",
        name = symbol,
        meta = biotype,
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
