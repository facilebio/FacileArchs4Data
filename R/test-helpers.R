#' Return an example sample data frame from a archs4 dataset
#' @export
example_archs4_samples <- function(
    species = c("mouse", "human"),
    afds = NULL,
    ...) {
  if (is.null(afds)) {
    species <- match.arg(species, c("mouse", "human"))
    stopifnot(species == "mouse")
    afds <- FacileArchs4DataSet(species)
  } else {
    species <- organism(afds)
  }
  expect_class(afds, "FacileArchs4DataSet")
  stopifnot(species == "mouse")
  if (species == "mouse") {
    out <- samples(afds) |> filter(dataset == "GSE135539")
    descr <- sub("Sample \\d+: ", "", out$title)
    dparts <- strsplit(descr, " ")
    out <- out |>
      dplyr::transmute(
        h5idx,
        dataset,
        sample_id,
        description = descr,
        genotype = ifelse(grepl("wt", descr, ignore.case = TRUE), "WT", "SOD1"),
        treatment = dplyr::case_when(
          grepl("untreated", descr, ignore.case = TRUE) ~ "untreated",
          grepl("control", descr, ignore.case = TRUE) ~ "untreated",
          grepl("shrna", descr, ignore.case = TRUE) ~ "shRNA",
          grepl("sham", descr, ignore.case = TRUE) ~ "sham",
          .default = "wtf"),
        region = sapply(dparts, tail, 1),
        subject_id = sapply(dparts, \(x) paste0("m", x[length(x) - 1])),
        batch = ifelse(grepl("Batch 1", characteristics_ch1), "b1", "b2")
      ) |>
      dplyr::mutate(
        group = paste(region, genotype, treatment, sep = "__"),
        genotype = factor(genotype) |> relevel("WT"),
        treatment = factor(treatment, c("untreated", "sham", "shRNA")),
        .after = "sample_id") |>
      dplyr::arrange(region, genotype, treatment) |>
      dplyr::mutate(group = factor(group, unique(group)))
  } else {
    stop("No human curation yet")
  }
  out
}
