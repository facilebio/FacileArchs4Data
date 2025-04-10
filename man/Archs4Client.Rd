% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Archs4Client.R
\name{Archs4Client}
\alias{Archs4Client}
\title{Thin wrapper to the hdf5 files of an Archs4 dataset.}
\description{
Thin wrapper to the hdf5 files of an Archs4 dataset.

Thin wrapper to the hdf5 files of an Archs4 dataset.
}
\details{
Create a wrapper to the hdf5 files that are directly downloaded from the
archs4 website.

To make your life easy, store the "equally versioned" mouse and human
datasets in a single directory with the value of \code{species} in the filename.
You can put new/old h5 files in different directories, then load the dataset
by specifying simply the species. See \strong{Examples} for examples.
}
\examples{

## ------------------------------------------------
## Method `Archs4Client$new`
## ------------------------------------------------

\dontrun{
options(archs4.data_dir = "/path/to/directory/with/hdf5-count-files")
a4m <- Archs4Client$new("mouse")
a4h <- Archs4Client$new("/path/to/human_gene_v2.2.h5")
}
}
\seealso{
https://github.com/MaayanLab/archs4py/blob/main/archs4py/data.py#L27
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{path}}{the filepath to the species-specific hdf5 data object}

\item{\code{samples}}{the facile_frame of the samples in this object}

\item{\code{features}}{the tibble of features (genes) in this dataset.}

\item{\code{species}}{is this thing mouse or human?}

\item{\code{remove_sc}}{default value for the \code{remove_sc} parameters}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Archs4Client-new}{\code{Archs4Client$new()}}
\item \href{#method-Archs4Client-search}{\code{Archs4Client$search()}}
\item \href{#method-Archs4Client-studies}{\code{Archs4Client$studies()}}
\item \href{#method-Archs4Client-clone}{\code{Archs4Client$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Archs4Client-new"></a>}}
\if{latex}{\out{\hypertarget{method-Archs4Client-new}{}}}
\subsection{Method \code{new()}}{
Create a new Archs4Client object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Archs4Client$new(
  species,
  ...,
  path = NULL,
  directory = getOption("archs4.data_dir", NULL),
  remove_sc = TRUE,
  sample_columns = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{species}}{\code{"mouse"} or \code{"human"}}

\item{\code{...}}{stuff}

\item{\code{path}}{the full path to the h5 file to load}

\item{\code{directory}}{the parent directory that stores the species-specific
expression datasets. Defaults to \code{getOption("archs4.data_dir")}. To make
your life easier, you can set this option in your \verb{~/.Rprofile}}

\item{\code{remove_sc}}{default value for the \code{remove_sc} parameters}

\item{\code{sample_columns}}{wut?}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new Archs4Client
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
options(archs4.data_dir = "/path/to/directory/with/hdf5-count-files")
a4m <- Archs4Client$new("mouse")
a4h <- Archs4Client$new("/path/to/human_gene_v2.2.h5")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Archs4Client-search"></a>}}
\if{latex}{\out{\hypertarget{method-Archs4Client-search}{}}}
\subsection{Method \code{search()}}{
Find samples that match a regular expression in any of the descriptive
columns.

This function is fashioned after the archs4py.meta() function.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Archs4Client$search(
  regex,
  meta_fields = NULL,
  scrub_search = FALSE,
  ignore.case = TRUE,
  remove_sc = self$remove_sc,
  threshold_sc = 0.5
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{regex}}{the regular expression to use for search}

\item{\code{meta_fields}}{which columns in the dataset to search against, defaults
to a set of columns defined in the archs4py package}

\item{\code{scrub_search}}{remove whitespace,dashes,etc. from values in
samples columns before running the regex? Setting this to \code{TRUE} slows
down the search appreciable (~ 3x - 4x). Default: \code{FALSE}}

\item{\code{ignore.case}}{Should we ignore case when using \code{regex}?
Default: \code{TRUE}}

\item{\code{remove_sc}}{Wether to remove single cell samples from search?
Default: \code{TRUE}}

\item{\code{threshold_sc}}{the probability threshold to indicate wheter sample
is from a singlecell dataset, default: 0.5}

\item{\code{...}}{stuff and things}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the samples facile_frame that matches the search
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Archs4Client-studies"></a>}}
\if{latex}{\out{\hypertarget{method-Archs4Client-studies}{}}}
\subsection{Method \code{studies()}}{
Return a table of studies (series_id) in here, optionally dropping
those that look singlecell-like
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Archs4Client$studies(
  sample_search = NULL,
  meta_fields = .default_meta_fields(),
  scrub_search = FALSE,
  ignore.case = TRUE,
  threshold_sc = 0.5,
  threshold_sc_likely = 0.2
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sample_search}}{If a search string is provided, this will trigger
a sample search via the \verb{$search()} method. The study table returned
will be augmented with a \code{search_hits} column that tells you how many
samples within the study were a hit for your search.}

\item{\code{scrub_search, meta_fields, ignore.case}}{parameters passed to \code{search} function.}

\item{\code{threshold_sc}}{the probability threshold to indicate wheter sample
is from a singlecell dataset, default: 0.5}

\item{\code{threshold_sc_likely}}{the proportion of samples in a study that
must look like singlecell data in order to classify the study as
singlecell. default: 0.2}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a tibble with 1 row per geo series (GSE) / study with the
following columns:
\itemize{
\item \code{nsamples}: number of samples in the study
\item \code{nsc}: number of samples that look like single cell, given the
\code{threshold_sc} parameter
\item \code{sc_likely}: if the number of samples in the study that look like
they are singlecell exceeds \code{threshold_sc_likely}, this is set to
\code{TRUE}
\item \code{search_hits}: The number of samples in the series that were a hit
for the \code{sample_search} query. If no \code{sample_search} was provided,
this is just \code{0}.
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Archs4Client-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Archs4Client-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Archs4Client$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
