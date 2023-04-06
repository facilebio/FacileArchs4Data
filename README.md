
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileArchs4DataSet

This package implements the Facile API over the ARCHS4 (v2) datasets,
and provides some fundamental code required to build the children
dataset repositories for the human and mouse [HDF5
files](https://maayanlab.cloud/archs4/download.html) files provided by
the ARCHS4 project.

These [ARCHS4 dataset](https://maayanlab.cloud/archs4/index.html)
repositories can then be loaded into a `FacileArchs4DataSet` object,
which is a subclass of a `FacileDataSet` and subsequently used by the
[facilebio](https://facile.bio/) ecosystem.

## Installation

To build the species-specifc dataset packages, users will need to
download the (v2) gene level HDF5 files for the species of interest
(mouse or human) from the [ARCHS4
project](https://maayanlab.cloud/archs4/download.html) and run the
`build` command:

``` r
# pak::pkg_install("facilebio/FacileArchs4DataSet")
FacileArchs4DataSet::build(
  h5 = "/path/to/archs4_gene_human_v2.1.2.h5",
  outdir = "/path/to/repository/dir",
  species = "human")
```

This will build the `FacileArchs4HumanDataSet` package in the
`/path/to/package/parentdir/FacileArchs4HumanDataSet` directory.

## Data Exploration

Once the data packge is successfully built, you can use it for
exploratory analysis, the facile way:

``` r
library(FacileAnalysis)

a4h <- FacileArchs4DataSet::load("/path/to/repository/dir/human")
a4h |>
  filter_samples(dataset == "GSE89189") |>
  fdgeGadget()
```

Refer to the
[FacileAnalysis](https://facilebio.github.io/FacileAnalysis/)
documentation for details on how to easily perform diferential
expression analysis, GSEA, etc.

## License

In order for this package to be useful, the end user must download the
HDF5 files provifed by the [ARCHS4
project](https://maayanlab.cloud/archs4/download.html). Use of the
ARCHS4 data is subject to its own [terms of
use](https://maayanlab.cloud/archs4/help.html), and it is the end-user’s
responsibility to ensure that they are using it in a compliant manner.

The code in this package that facilitates the query and retrieval of the
data provided is released under [The Apache 2
License](https://www.apache.org/licenses/LICENSE-2.0).

## Acknowledgements

Thanks toAlexander Lachmann and The Ma’ayan Laboratory for the
development andcontinued updates to the [ARCHS4
project](https://maayanlab.cloud/archs4/index.html).
