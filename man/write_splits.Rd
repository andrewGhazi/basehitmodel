% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_proteins_separately.R
\name{write_splits}
\alias{write_splits}
\title{Write split data}
\usage{
write_splits(filtered_data, split_data_dir, verbose = TRUE)
}
\description{
Write out the data for each protein separately. It's written out
  to files named by index because sometimes proteins have weird characters or
  slashes in their names.
}
