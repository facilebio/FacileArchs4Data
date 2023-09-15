#' @noRd
#' @export
biocbox.archs4_client_facile_frame <- function(
    x, class = "list", ...) {
  a4c <- assert_class(fds(x), "Archs4Client")
  class <- assert_choice(class, c("list", "DGEList", "SummarizedExperiment"))

  features <- a4c$features
  counts <- rhdf5::h5read(
    a4c$path,
    "data/expression",
    index = list(x$h5idx, features$h5idx)) |>
    t()
  rownames(counts) <- features$ensembl_gene_id
  colnames(counts) <- x$sample

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
    samples = x |>
      as.data.frame() |>
      select(-h5idx) |>
      rename(sample_id = sample, dataset = series_id))
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
