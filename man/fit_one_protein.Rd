% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_proteins_separately.R
\name{fit_one_protein}
\alias{fit_one_protein}
\title{Fit the model to one protein}
\usage{
fit_one_protein(
  protein,
  split_data_dir,
  save_fits,
  save_summaries,
  proteins,
  ixn_prior_width,
  iter_sampling,
  iter_warmup = 1000,
  algorithm,
  out_dir
)
}
\arguments{
\item{protein}{the protein}

\item{p1_eta}{the subset of bh_eta for the protein in question,}

\item{stan_model}{protein model to use}
}
\description{
Fit the model to one protein
}
