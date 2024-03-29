% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_model_inputs}
\alias{plot_model_inputs}
\title{Plot barcode level model inputs}
\usage{
plot_model_inputs(
  bh_input,
  barcodes = NULL,
  proteins = NULL,
  strains = NULL,
  force = FALSE,
  log10_counts = TRUE
)
}
\arguments{
\item{bh_input}{a data frame of filtered model inputs}

\item{barcodes}{a character vector of barcodes to filter to}

\item{proteins}{a character vector of proteins to filter to}

\item{strains}{a character vector of strains to filter to}

\item{force}{if TRUE, override the check preventing gigantic plots}

\item{log10_counts}{if TRUE, log10 the values in both panels}
}
\description{
Plot barcode level model inputs
}
\details{
Grey cells in the top panel correspond to zeros in the original input (if log10_counts =
  TRUE). Empty cells (i.e. where you can see the underlying grid lines) are non-present in the
  input, likely due to being filtered out.
}
\examples{
\dontrun{
bh_input = data.table::fread("~/Desktop/tmp/cache/bh_input.tsv.gz")
plot_model_inputs(bh_input,
                  strains = c("AB1", "AB10", "AB12"),
                  proteins = c("LSAMP", "THSD1_Epitope-1", "LRTM1"))
}
}
