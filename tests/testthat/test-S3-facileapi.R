if (!exists("afds")) afds <- FacileArchs4DataSet("mouse")
if (!exists("asamples"))asamples <- example_archs4_samples(afds = afds)

test_that("biocbox returns matrix in prespecified order", {
  bb.all <- biocbox(asamples, class = "list")

  # scramle gene order, and make sure it comes back in the same way
  ff <- features(afds) |> dplyr::sample_n(10)
  bb.ff <- biocbox(asamples, features = ff, class = "list")
  expect_equal(ff$feature_id, rownames(bb.ff$assay_data))
  expect_equal(ff$feature_id, rownames(bb.ff$features))

  expect_equal(
    bb.ff$assay_data,
    bb.all$assay_data[ff$feature_id,],
    check.attributes = FALSE
  )

  # scramble sample order
  ss <- asamples |> dplyr::sample_n(10)
  bb.ss <- biocbox(ss, features = ff, class = "list")
  expect_equal(
    bb.ss$assay_data,
    bb.ff$assay_data[,ss$sample_id],
    check.attributes = FALSE
  )
})

test_that("(fetch|with)_assay_data works", {
  ff <- features(afds) |> sample_n(5)
  exprs <- fetch_assay_data(asamples, ff)
  nexprs <- fetch_assay_data(asamples, ff, normalized = TRUE)
  scores <- fetch_assay_data(asamples, ff, normalized = TRUE, aggregate = TRUE)
})

test_that("fetch_sample_covariates works on a skinny table", {
  xs <- select(asamples, dataset, sample_id)
  xslong <- fetch_sample_covariates(xs, "title")
  expect_character(xslong$variable, fixed = "title")
  expect_character(xslong$value)
  expect_character(xslong$class, fixed = "categorical")

  xswide <- with_sample_covariates(xs, c("title", "alignedreads"))
  expect_character(xswide$title)
  expect_numeric(xswide$alignedreads)
})

test_that("fetch_sample_covariates works on skinny table from duckdb archs4", {
  ddb <- FacileArchs4DataSet("mouse", .duckdb = TRUE)
  sid <- asamples$dataset[1L]
  xs <- samples(ddb) |>
    filter(.data$dataset == sid) |>
    select(dataset, sample_id)

  xslong <- fetch_sample_covariates(xs, "title")
  expect_character(xslong$variable, fixed = "title")
  expect_character(xslong$value)
  expect_character(xslong$class, fixed = "categorical")

  xswide <- with_sample_covariates(xs, c("title", "alignedreads"))
  expect_character(xswide$title)
  expect_numeric(xswide$alignedreads)
})
