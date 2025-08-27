devtools::load_all(".")

if (!exists("afds")) afds <- FacileArchs4DataSet("mouse")
if (!exists("asamples"))asamples <- example_archs4_samples(afds = afds)

library(FacileAnalysis)
flm <-  flm_def(asamples, "batch")
xdge <- fdge(flm, method = "voom")

# debug(FacileAnalysis:::flm_def.facile_frame)
# debug(fetch_sample_covariates.FacileArchs4DataSet)


FacileAnalysisShine::fdgeGadget(asamples)


samples.all <- FacileReviRData::load_samples("rtx117-vwm-xtalpi") |>
  dplyr::mutate(stupid = genotype)

FacileAnalysisShine::fdgeGadget(samples.all)
