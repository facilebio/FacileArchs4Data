# Implementing same functions required in FacileBiocData to make some magic

#' @noRd
#' @export
biocbox.archs4_facile_frame <- function(
    x, class = "list", assay_name = "counts", features = NULL,
    sample_covariates = NULL, with_h5idx = FALSE, ...) {
  assert_class(x, "archs4_facile_frame")
  assert_integerish(x$h5idx)
  afds <- assert_class(FacileData::fds(x), "FacileArchs4DataSet")
  class <- assert_choice(class, c("list", "DGEList", "SummarizedExperiment"))

  if (is.null(features)) {
    features <- features(afds)
  }
  assert_tibble(features)
  assert_subset(c("h5idx", "feature_id", "name", "meta"),  names(features))
  assert_integerish(features$h5idx)
  assert_character(features$feature_id)
  stopifnot("duplicated features" = !any(duplicated(features$feature_id)))

  counts <- rhdf5::h5read(
    m4$h5,
    "data/expression",
    index = list(x$h5idx, features$h5idx)
  ) |>
    t()
  rownames(counts) <- features$feature_id

  genes <- as.data.frame(features) |>
    dplyr::mutate(feature_type = "esngid", source = "ensembl_v107")
  rownames(genes) <- rownames(counts)

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

  samples <- as.data.frame(x)
  rownames(samples) <- colnames(counts)

  out <- list(counts = counts, features = genes, samples = samples)

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

#' @noRd
#' @export
samples.FacileArchs4DataSet <- function(x, ..., drop_covariates = FALSE) {
  assert_class(x, "FacileArchs4DataSet")
  assert_flag(drop_covariates)
  out <- x$samples
  if (drop_covariates) {
    out <- dplyr::select(out, dataset, sample_id)
  }
  out
}

#' @noRd
#' @export
features.FacileArchs4DataSet <- function(x, ...) {
  assert_class(x, "FacileArchs4DataSet")
  x$features
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
      h5idx,
      assay = .env$assay_name,
      libsize = .data$alignedreads,
      normfactor = 1,
      singlecellprobability
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
#' @examples
#' m4 <- FacileArchs4DataSet()
#' xf <- features(m4) |> filter(name %in% c("Sesn2", "Fgf21"))
#' xs <- samples(m4) |> filter(grepl("GSE249764", dataset))
fetch_assay_data.FacileArchs4DataSet <- function(
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
    prior.count = 0.25,
    subset.threshold = 700,
    aggregate = FALSE,
    aggregate.by = "ewm",
    verbose = FALSE
) {
  assert_flag(as.matrix)
  assert_flag(normalized)
  if (is.null(assay_name)) assay_name <- default_assay(x)
  assert_choice(assay_name, assay_names(x))

  aggregate.by <- match.arg(tolower(aggregate.by), c("ewm", "zscore"))
  ainfo <- assay_info(x, assay_name)

  samples.all <- samples(x)
  if (is.null(samples)) {
    samples <- samples.all
  } else {
    assert_sample_subset(samples)
    minfo <- dplyr::select(samples.all, dataset, sample_id, h5idx)
    samples <- samples |>
      dplyr::inner_join(
        minfo,
        by = c("dataset", "sample_id"),
        suffix = c(".insample", "")
      )
  }
  assert_data_frame(samples)
  samples <- samples |>
    dplyr::mutate(samid = paste(dataset, sample_id, sep = "__"))

  assert_integerish(samples$h5idx)

  features.all <- features(x, assay_name)
  if (is.null(features)) {
    features <- features.all
  } else {
    if (is.data.frame(features)) {
      features <- features[["feature_id"]]
    }
    if (is.factor(features)) features <- as.character(features)
    if (is.character(features)) {
      features <- filter(features.all, .data$feature_id %in% .env$features)
    }
    assert_data_frame(features)
    if (!"assay" %in% colnames(features) || !is.character(features$assay)) {
      features <- collect(features, n = Inf)
      features[["assay"]] <- assay_name
    }
    assert_assay_feature_descriptor(features)
  }
  assert_data_frame(features)
  assert_integerish(features$h5idx)
  features <- dplyr::distinct(features, .data$feature_id, .keep_all = TRUE)
  features[["assay_type"]] <- ainfo[["assay_type"]]

  adat <- rhdf5::h5read(
    x$h5,
    "data/expression",
    index = list(samples$h5idx, features$h5idx)
  ) |>
    t()
  rownames(adat) <- features$feature_id
  colnames(adat) <- samples$samid

  if (!is.null(batch)) {
    if (!isTRUE(normalized)) {
      warning("`batch` parameter specified, setting `normalized` to `TRUE`")
    }
    normalized <- TRUE
  }

  if (normalized) {
    # Adds sample-level assay data appropriate for whatever the assay is
    asinfo <- assay_sample_info(x, assay_name, samples)
    # If `samples` were passed in with any assay-level covariates, let those
    # override what is in the database
    custom.cols <- intersect(colnames(samples), colnames(asinfo))
    custom.cols <- setdiff(custom.cols, c("dataset", "sample_id"))
    if (length(custom.cols)) {
      asinfo <- asinfo[, !colnames(asinfo) %in% custom.cols]
    }
    if (length(setdiff(colnames(asinfo), c("dataset", "sample_id")))) {
      samples <- left_join(samples, asinfo, by = c("dataset", "sample_id"))
    }
    adat <- normalize_assay_data(adat, features, samples, batch = batch,
                                 main = main, verbose = verbose, ...)
  }

  aggregated <- NULL

  if (isTRUE(aggregate)) {
    if (aggregate.by == "ewm") {
      aggregated <- sparrow::eigenWeightedMean(adat, ...)
    } else if (aggregate.by == "zscore") {
      aggregated <- sparrow::zScore(adat, ...)
    } else {
      stop("Unknown aggregation method: ", aggregate.by)
    }
    adat <- matrix(
      aggregated$score,
      nrow = 1L,
      dimnames = list("score", names(aggregated$score)))
  }

  if (!as.matrix) {
    atype <- ainfo[["assay_type"]]
    ftype <- ainfo[["feature_type"]]
    vals <- FacileData:::.melt.assay.matrix(
      adat,
      assay_name,
      atype,
      ftype,
      features
    )
    vals <- as_tibble(vals)
    if (isTRUE(aggregate)) {
      vals <- mutate(
        vals,
        feature_type = "aggregated",
        feature_id = "aggregated",
        feature_name = "aggregated"
      )
    }
    vals <- samples |>
      left_join(
        vals,
        by = c("dataset", "sample_id"),
        suffix = c("", ".dropme..")) |>
      select(-ends_with(".dropme.."))
  } else {
    vals <- adat
  }

  set_fds(vals, x)
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
