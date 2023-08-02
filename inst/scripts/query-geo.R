# datasets of interest related to Legouis et al. 2020
# library(rhdf5)
# library(dplyr)

devtools::load_all(".")
leg.ids <- c(
  "GSE151167", # mouse RNAseq, 5 samples (2020)
  "GSE52004",  # mouse affy chip, 57 samples (2014)
  "GSE126805", # human rnaseq, 163 samples
  NULL)

h5 <- list(
  human = "~/workspace/data/archs4/v2.2/human_gene_v2.2.h5",
  mouse = "~/workspace/data/archs4/v2.2/mouse_gene_v2.2.h5")


a4hs <- Archs4Client$new(h5$human)
a4mm <- Archs4Client$new(h5$mouse)

genes <- .hdf5_group_load_table(h5$human, "meta/genes")
samples <- .hdf5_group_load_table(h5$human, "meta/samples")


