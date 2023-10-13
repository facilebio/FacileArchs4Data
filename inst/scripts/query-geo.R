# datasets of interest related to Legouis et al. 2020
# library(rhdf5)
library(dplyr)

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


GSE126805 <- filter(a4hs$samples, series_id == "GSE126805")
GSE151167 <- filter(a4mm$samples, series_id == "GSE151167")
GSE149739 <- filter(a4mm$samples, series_id == "GSE149739")

se <- biocbox(GSE149739, "SummarizedExperiment")
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# data from Humphreys lab, should have 59 samplesn across mouse and human
GSE199437 <- list(
  mouse = filter(a4mm$samples, series_id == "GSE199437"),  # 46
  human = filter(a4hs$samples, series_id == "GSE199437"))  # 9
sapply(GSE199437, nrow) |> sum() # 55

