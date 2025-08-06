if (!exists("afds")) afds <- FacileArchs4DataSet("mouse")
if (!exists("asamples")) {
  asamples <- example_archs4_samples(afds = afds)
}

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
