# Implementing same functions required in FacileBiocData to make some magic

#' @noRd
#' @export
name.FacileArchs4DataSet <- function(x, ...) {
  org <- gsub(" +", "", tools::toTitleCase(organism(x)))
  sprintf("FacileArchs4%sDataSet", tools::toTitleCase(organism(x)))
}

#' @noRd
#' @export
organism.FacileArchs4DataSet <- function(x, ...) {
  x$species
}

#' @noRd
#' @export
biocbox.archs4_facile_frame <- function(
  x,
  class = "list",
  assay_name = "counts",
  features = NULL,
  sample_covariates = NULL,
  normalized = FALSE,
  with_h5idx = FALSE,
  ...
) {
  x <- assert_class(x, "archs4_facile_frame") |> dplyr::collect()
  assert_character(x$dataset)
  assert_character(x$sample_id)
  assert_integerish(x$h5idx)
  afds <- assert_class(FacileData::fds(x), "FacileArchs4DataSet")
  class <- assert_choice(class, c("list", "DGEList", "SummarizedExperiment"))

  if (is.null(features)) {
    features <- features(afds)
  } else {
    features <- .extract_features(features, afds)
  }
  assert_tibble(features)
  assert_subset(c("h5idx", "feature_id", "name", "meta"), names(features))
  assert_integerish(features$h5idx)
  assert_character(features$feature_id)
  stopifnot("duplicated features" = !any(duplicated(features$feature_id)))

  # counts <- rhdf5::h5read(
  #   m4$h5,
  #   "data/expression",
  #   index = list(x$h5idx, features$h5idx)
  # ) |>
  #   t()
  # rownames(counts) <- features$feature_id

  A <- fetch_assay_data(
    afds,
    features,
    x,
    assay_name = assay_name,
    normalized = normalized,
    as.matrix = TRUE,
    ...
  )
  stopifnot(
    "feature order check" = all.equal(rownames(A), features$feature_id),
    "sample order check" = all.equal(
      paste(x$dataset, x$sample_id, sep = "__"),
      colnames(A)
    )
  )

  genes <- as.data.frame(features) |>
    dplyr::mutate(
      symbol = name,
      feature_type = "esngid",
      source = "ensembl_v107"
    )
  rownames(genes) <- genes$feature_id

  # if (!any(c("sample_id", "sample") %in% colnames(x))) {
  #   warning(
  #     "neither `sample_id` or `sample` column found in sample frame, ",
  #     "creating pseudo id's"
  #   )
  #   x$sample_id <- sprintf("sample_%d", seq_len(ncol(A)))
  # }
  # if ("sample_id" %in% colnames(x)) {
  #   colnames(A) <- x$sample_id
  # } else {
  #   colnames(A) <- x$sample
  # }

  # Let's work with the sample frame we have, not the one we think we have, ie.
  # we may have manipulated / cleaned this up before passing it in here.
  # if (!"dataset" %in% names(x) && "series_id" %in% names(x)) {
  #   x <- dplyr::rename(x, dataset = series_id)
  # }
  if (!with_h5idx) {
    x$h5idx <- NULL
  }

  samples <- as.data.frame(x)
  rownames(samples) <- colnames(A)

  out <- list(assay_data = A, features = genes, samples = samples)

  if (class == "DGEList") {
    reqpkg("edgeR")
    out <- edgeR::DGEList(
      counts = out$assay_data,
      genes = out$features,
      samples = out$samples
    ) |>
      edgeR::calcNormFactors()
  } else if (class == "SummarizedExperiment") {
    reqpkg("SummarizedExperiment")
    out <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = out$A),
      rowData = out$features,
      colData = out$samples
    )
  }

  out <- FacileData::set_fds(out, afds)
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
features.FacileArchs4DataSet <- function(
  x,
  assay_name = "counts",
  feature_type = "ensgid",
  feature_ids = NULL,
  ...
) {
  assert_class(x, "FacileArchs4DataSet")
  stopifnot(
    "only counts assay" = identical(assay_name, "counts"),
    "only ensembl gens" = identical(feature_type, "ensgid")
  )
  out <- x$features
  if (!is.null(feature_ids)) {
    out <- .extract_features(feature_ids, x)
  }
  out
}

# assay data methods -----------------------------------------------------------

#' @noRd
#' @export
#' @examples
#' mfds <- FacileArchs4DataSet("mouse")
#' xs <- samples(mfds) |> filter(dataset == "GSE116863")
#' asi <- assay_sample_info(mfds, samples = xs)
#'
assay_sample_info.FacileArchs4DataSet <- function(
  x,
  assay_name = "counts",
  samples = NULL,
  ...,
  verbose = FALSE,
  conflict_drop = c("samples", "fds"),
  .developer = getOption("fbioc.developer", FALSE)
) {
  conflict_drop <- match.arg(conflict_drop)
  conflict_keep <- setdiff(c("samples", "fds"), conflict_drop)
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
        suffix = c(".samples", ".fds")
      )

    conflicted <- colnames(out)[grepl("\\.(samples|fds)$", colnames(out))]
    if (length(conflicted) > 0L) {
      if (verbose) {
        warning(
          "conflicting assay column information in samples and fds, will drop",
          "values from `", conflict_drop, "`: ",
          paste(conflicted, collapse = ",")
        )
      }
      out <- dplyr::select(out, -dplyr::ends_with(conflict_drop))
      colnames(out) <- sub(paste0("\\.", conflict_keep, "$"), "", colnames(out))
      noinfo <- dplyr::filter(out, is.na(libsize))
      if (nrow(noinfo) > 0L) {
        warning(
          "These samples were not found in archs4 dataset: ",
          paste(noinfo$sample_id, collapse = ",")
        )
      }
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
    ~assay,
    ~assay_type,
    ~feature_type,
    ~description,
    ~nfeatures,
    ~storage_mode,
    "counts",
    "rnaseq",
    "ensgid",
    "counts from archs4",
    nrow(x$features),
    "integer"
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
    features <- .extract_features(features, x)
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
    adat <- normalize_assay_data(
      adat,
      features,
      samples,
      batch = batch,
      main = main,
      verbose = verbose,
      ...
    )
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
      dimnames = list("score", names(aggregated$score))
    )
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
    vals <- dplyr::as_tibble(vals)
    if (isTRUE(aggregate)) {
      vals <- dplyr::mutate(
        vals,
        feature_type = "aggregated",
        feature_id = "aggregated",
        feature_name = "aggregated"
      )
    }
    vals <- samples |>
      dplyr::left_join(
        vals,
        by = c("dataset", "sample_id"),
        suffix = c("", ".dropme..")
      ) |>
      dplyr::select(-dplyr::ends_with(".dropme.."))
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

# sample covariate methods -----------------------------------------------------
# by only defining fetch_sample_covariates.FacileArchs4DataStore, we don't
# have to define the with_* functions

#' @noRd
#' @export
fetch_sample_covariates.FacileArchs4DataSet <- function(
  x,
  covariates = NULL,
  samples = NULL,
  custom_key = Sys.getenv("USER"),
  with_source = FALSE,
  ...
) {
  sinfo <- samples(x)
  if (is.null(samples)) {
    stop("You don't want to call fetch_sample_covariates on all of archs4")
  } else {
    samples <- assert_sample_subset(dplyr::collect(samples))
    sinfo <- sinfo |>
      dplyr::filter(
        .data$dataset %in% samples$dataset & .data$sample_id %in% samples$sample_id
      ) |>
      dplyr::collect()
    sinfo <- dplyr::semi_join(sinfo, samples, by = c("dataset", "sample_id"))
    if (nrow(sinfo) == 0L) {
      return(.empty_sample_covariates(x))
    }
  }
  all.covs <- setdiff(colnames(sinfo), c("dataset", "sample_id"))

  if (is.null(covariates)) {
    covariates <- all.covs
  } else {
    assert_character(covariates)
    covariates <- setdiff(covariates, c("dataset", "sample_id"))
    # browser()
    bad.covs <- setdiff(covariates, all.covs)
    if (length(bad.covs)) {
      warning("Unknown sample covariates: ", paste(bad.covs, collapse = ","))
      covariates <- intersect(covariates, all.covs)
    }
    sinfo <- dplyr::select(sinfo, "dataset", "sample_id", {{ covariates }})
  }

  out <- FacileData::as.EAVtable(sinfo)
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @noRd
#' @export
covariate_definitions.FacileArchs4DataSet <- function(x, as.list = TRUE, ...) {
  # eav <- as.EAVtable(pdata(x))
  # out <- attr(eav, "covariate_def")
  # out <- ifacile(x)[["covariate_def"]]

  # sample table is huge, and it's not so informative, so just recalculatre
  # eav_metadata based on first few columns
  out <- samples(x) |>
    dplyr::collect(n = 5) |>
    eav_metadata_create()
  if (!as.list) {
    out <- lapply(names(out), function(name) {
      i <- out[[name]]
      lvls <- unique(i$levels)
      is.factor <- !is.null(lvls)
      lbl <- if (is.null(i$label)) name else i$label
      tibble(
        variable = name,
        type = i$type,
        class = i$class,
        label = i$label,
        is_factor = is.factor,
        levels = list(lvls),
        description = i$description
      )
    })
    out <- bind_rows(out)
  }
  class(out) <- c('CovariateDefinitions', class(out))
  set_fds(out, x)
}

#' @noRd
#' @export
fetch_custom_sample_covariates.FacileArchs4DataSet <- function(
  x,
  covariates = NULL,
  samples = NULL,
  custom_key = Sys.getenv("USER"),
  file.prefix = "facile",
  ...
) {
  .empty_sample_covariates(x)
}

#' @noRd
.empty_sample_covariates <- function(x = NULL) {
  out <- tibble(
    dataset = character(),
    sample_id = character(),
    variable = character(),
    value = character(),
    class = character(),
    type = character(),
    date_entered = integer()
  )
  if (is(x, "FacileArchs4DataSet")) {
    out <- as_facile_frame(out, x, .valid_sample_check = FALSE)
  }
  out
}

# misc utilities ---------------------------------------------------------------

#' Extract the feature descriptor from the archs4 dataset
.extract_features <- function(x, afds) {
  assert_class(afds, "FacileArchs4DataSet")
  features.all <- afds$features
  if (is.data.frame(x)) {
    x <- x[["feature_id"]]
  }
  if (is.factor(x)) x <- as.character(x)
  assert_character(x)
  x <- dplyr::tibble(feature_id = x) |>
    dplyr::left_join(features.all, by = "feature_id")
  missed <- dplyr::filter(x, is.na(h5idx))
  if (nrow(missed) > 0) {
    stop("Missing arcsh4 features: ", paste(missed$feature_id, collapse = ","))
  }
  assert_data_frame(x)
  if (!"assay" %in% colnames(x) || !is.character(x$assay)) {
    if (!length(assay_names(afds)) == 1) {
      stop("This needs updating if we want to support multi-assay archs4")
    }
    x[["assay"]] <- "counts"
  }
  (assert_assay_feature_descriptor(x))
}
