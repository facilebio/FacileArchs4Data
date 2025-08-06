# Implementing same functions required in FacileBiocData to make some magic

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

#' @noRd
#' @export
#' @examples
#' mfds <- FacileArchs4DataSet("mouse")
#' xs <- samples(mfds) |> filter(dataset == "GSE116863")
#' asi <- assay_sample_info(mfds, samples = xs)
#'
assay_sample_info.FacileArchs4DataSet  <- function(
    x, assay_name = "counts", samples = NULL, ...,
    .developer = getOption("fbioc.developer", FALSE)) {
  assert_choice(assay_name, assay_names(x))
  if (assay_name != "counts") {
    warning(
      "We haven't figured out where to get libisze and normfactor for ",
      "anything but counts assay"
    )
  }

  asi.all <- samples(x, drop_covariates = FALSE) |>
    dplyr::transmute(
      dataset,
      sample_id,
      assay = .env$assay_name,
      libsize = .data$alignedreads,
      normfactor = 1,
      h5idx
    )

  if (is.null(samples)) {
    out <- asi.all
  } else {
    out <- samples |>
      assert_sample_subset() |>
      left_join(
        asi.all,
        by = c("dataset", "sample_id"),
        suffix = c(".samples", "")
      )
    conflicted <- colnames(out)[endsWith(colnames(out), ".samples")]
    if (length(conflicted) > 1L) {
      waring(
        "assay information colums already found in samples frame: ",
        paste(conflicted, collapse = ","))
    }
    noinfo <- dplyr::filter(out, is.na(libsize))
    if (nrow(noinfo) > 0L) {
      warning(
        "These samples were not found in archs4 dataset: ",
        paste(noinfo$sample_id, collapse = ","))
    }
  }
  out
}


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
default_assay.FacileArchs4DataSet <- function(x, ...) {
  "counts"
}

#' @noRd
#' @export
assay_names.FacileArchs4DataSet <- function(x, default_first = TRUE, ...) {
  "counts"
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
