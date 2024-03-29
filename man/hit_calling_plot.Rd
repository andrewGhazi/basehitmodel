% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{hit_calling_plot}
\alias{hit_calling_plot}
\title{Make a hit calling plot}
\usage{
hit_calling_plot(
  protein_,
  strain_,
  bead_binding,
  concordances,
  fit_summary,
  weak_score_threshold = 0.5,
  strong_score_threshold = 1,
  weak_concordance = 0.75,
  strong_concordance = 0.95,
  tags = TRUE,
  intervals = c(0.95, 0.99)
)
}
\value{
a list of three plots
}
\description{
Show both the posterior intervals and placement in the
  concordance distribution
}
\details{
The concordance distribution shows the concordance values for all
  interactions that exclude 0 as a credible value by 95% interval
}
