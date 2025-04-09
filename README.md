
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileArchs4Data

This package is a work in progress and aims to implement the Facile API
over the ARCHS4 (v2) datasets. It provides some fundamental code
required to build the children dataset repositories for the human and
mouse [HDF5 files](https://maayanlab.cloud/archs4/download.html) files
provided by the ARCHS4 project.

These [ARCHS4 dataset](https://maayanlab.cloud/archs4/index.html)
repositories can then be loaded into a `FacileArchs4Data` object,
which is a subclass of a `FacileDataSet` and subsequently used by the
[facilebio](https://facile.bio/) ecosystem.

This package is currently written against
[v2.2](https://maayanlab.cloud/archs4/download.html) of the human and
mouse file format.

## State of Affairs

Currently you should use the `a4 <- Archs4Client$new("mouse")` (or
`"human"`) R6 constructor to get a thin wrapper around the data object.
Then, subset the `a4$samples` tibble to include the the samples you want
to analyze. Finally, convert the sample-frame it to a `DGEList` or
`SummarizedExperiment` to move on with your life – see the **Usage**
section for details.

To make your life easy, store the “equally versioned” mouse and human
datasets in a single directory (`$ARCHS4_DATADIR`) with the value of
`species` in the filename, and set a global R option (`archs4.data_dir`)
to `$ARCHS4_DATADIR`.

You can put new/old h5 files in different directories, then load the
dataset by specifying simply the species. See the **Installation**
section for more details.

## Usage

``` r
library(FacileArchs4Data)
# NOTE: Set this in your ~/.Rprofile file to make your life easier
options(archs4.data_dir = "~/workspace/data/archs4/v2.2")
a4h <- Archs4Client$new("human")
#> Warning in FUN(X[[i]], ...): NAs introduced by coercion to integer range
a4h$samples |> head()
#> # A tibble: 6 × 16
#>   h5idx series_id geo_accession sample     source_name_ch1      title           
#>   <int> <chr>     <chr>         <chr>      <chr>                <chr>           
#> 1     1 GSE29282  GSM1000981    GSM1000981 Human DLBCL cel line OCI-LY1_48hrs_m…
#> 2     2 GSE29282  GSM1000982    GSM1000982 Human DLBCL cel line OCI-LY1_48hrs_m…
#> 3     3 GSE29282  GSM1000983    GSM1000983 Human DLBCL cel line OCI-LY1_48hrs_m…
#> 4     4 GSE29282  GSM1000984    GSM1000984 Human DLBCL cel line OCI-LY1_48hrs_m…
#> 5     5 GSE29282  GSM1000985    GSM1000985 Human DLBCL cel line OCI-LY1_48hrs_m…
#> 6     6 GSE29282  GSM1000986    GSM1000986 Human DLBCL cel line OCI-LY1_48hrs_m…
#> # ℹ 10 more variables: singlecellprobability <dbl>, readsaligned <int>,
#> #   molecule_ch1 <chr>, library_source <chr>, library_selection <chr>,
#> #   platform_id <chr>, instrument_model <chr>, data_processing <chr>,
#> #   extract_protocol_ch1 <chr>, characteristics_ch1 <chr>
```

## Installation

Create a directory (`$ARCHS4_DATADIR`) on your local filesystem and
download the mouse and human gene expression count matrices from the
[ARCHS4 Download page](https://maayanlab.cloud/archs4/download.html).
Let’s assume `ARCHS4_DATA_DIR='~/data/archs4/v2.2'`

This directory should look like:

``` bash
$ ls -l $ARCHS4_DATA_DIR # ls -l ~/data/archs4/v2.2
-rw-r--r--@ 1 steve  staff    35G Jan 22 14:25 human_gene_v2.2.h5
-rw-r--r--@ 1 steve  staff    28G Feb  5 16:13 mouse_gene_v2.2.h5
```

In your ~/.Rprofile, set the `archs4.data_dir` option to
`$ARCHS4_DATADIR`:

``` bash
$ cat ~/.Rprofile
options(
  repos = c(cloud = "https://cloud.r-project.org"),
  archs4.data_dir = "~/data/archs4/v2.2")
```

Now install this package from GitHub.

``` r
remotes::install_github("facilebio/FacileArchs4Data")
```

## Usage

Find a human dataset of interest and convert it into a
`SummarizedExperiment`:

``` r
library(FacileArchs4Data)
a4h <- Archs4Client$new("human")
#> Warning in FUN(X[[i]], ...): NAs introduced by coercion to integer range
# Let's get a DGEList of some iPSC microglia data from GSE186301
y <- a4h$samples |> 
  dplyr::filter(series_id == "GSE186301") |> 
  biocbox(class = "DGEList")
y$samples |> 
  dplyr::select(dataset, sample_id, source_name_ch1) |> 
  dplyr::sample_n(5)
#>              dataset  sample_id source_name_ch1
#> GSM5643402 GSE186301 GSM5643402  iPSC-microglia
#> GSM5643403 GSE186301 GSM5643403  iPSC-microglia
#> GSM5643405 GSE186301 GSM5643405  iPSC-microglia
#> GSM5643407 GSE186301 GSM5643407  iPSC-microglia
#> GSM5643404 GSE186301 GSM5643404  iPSC-microglia
```

NOTE: The rest of these instructions was from v1 of this software. There
were some good ideas there, so leaving it for reference but they don’t
work.

## Installation

To build the species-specifc dataset packages, users will need to
download the (v2) gene level HDF5 files for the species of interest
(mouse or human) from the [ARCHS4
project](https://maayanlab.cloud/archs4/download.html) and run the
`build` command:

``` r
# pak::pkg_install("facilebio/FacileArchs4Data")
FacileArchs4Data::build(
  h5 = "/path/to/archs4_gene_human_v2.1.2.h5",
  outdir = "/path/to/repository/dir",
  species = "human")
```

This will create a `/path/to/repository/dir/human` directory that can be
used to instantiate a `FacileArcsh4DataSet` that can be used in the
facilebio ecosystem.

## Data Exploration

Once the data packge is successfully built, you can use it for
exploratory analysis, the facile way:

``` r
library(FacileAnalysis)

a4h <- FacileArchs4Data::load("/path/to/repository/dir/human")
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
